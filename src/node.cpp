/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "src/node.hpp"
#include "src/lattice.hpp"

class Lattice;

const int Node::dim  = 2;

const int Node::dir = 9;

const std::vector<double> 
Node::weights = {4./9., 
                1./9., 1./9., 1./9., 1./9., 
                1./36., 1./36., 1./36., 1./36.};

const std::vector<std::vector<double>> 
Node::coeff = {{0., 0.},
              {1., 0.}, {0., -1.}, {-1., 0.}, {0., 1.},
              {1., -1.}, {-1., -1.}, {-1., 1.}, {1., 1.}};
// 0 1 2 3 4 5 6 7 8 
const std::vector<int> 
Node::bb_indexes = {0, 3, 4, 1, 2, 7, 8, 5, 6};

Node::Node(NodeType node_type_, std::vector<unsigned int> coord_,
        double ux_, double uy_, double rho_)
  : node_type(node_type_),
    coord(coord_),
    ux(ux_),
    uy(uy_),
    rho(rho_),
    zou_he_type(ZouHeType::none),
{
  f_pre = std::make_unique<std::vector<double>>(dir, 0.);
  f_post = std::make_unique<std::vector<double>>(dir, 0.);
  f_adj = std::make_unique<std::vector<double>>(dir, 0.);
  bounce_back_dir.resize(dir, false);
  bounce_back_delta.resize(dir, 0.);
}

void
Node::init_equilibrium()
{
  // initialize f pre-collision to equilibrium distribution (only at step 0)
  // f_i = w_i * rho * (1  + 3 * (c_i.u) + 4.5 * (c_i.u)^2 - 1.5 * u^2)

  for(unsigned int i = 0; i<dir; ++i){
    (*f_pre)[i] = weights[i] * rho * (1. + 3.0 * (coeff[i][0] * ux + coeff[i][1] * uy) + 
              4.5 * (coeff[i][0] * ux + coeff[i][1] * uy) * (coeff[i][0] * ux + coeff[i][1] * uy) - 
              1.5 * (ux * ux + uy * uy));
  }
}

void
Node::collide(const Lattice &lattice)
{
  double tau = lattice.get_tau();
  double dt = lattice.get_dt();
  const double tauinv = 1./tau;   
  const double omtauinv = 1.0-tauinv;
  
  // Collision and streaming step: now compute and relax to equilibrium 
  // BGK collision operator
  // f_i(t+dt) = f_i(t) - dt/tau (f_i(t) - f_i^eq(t))

  double f_eq = 0.0;
  for(unsigned int i = 0; i<dir; ++i){
    f_eq =  weights[i] * rho * (1. + 3.0 * (coeff[i][0] * ux + coeff[i][1] * uy) + 
            4.5 * (coeff[i][0] * ux + coeff[i][1] * uy) * (coeff[i][0] * ux + coeff[i][1] * uy) - 
            1.5 * (ux * ux + uy * uy));
    (*f_post)[i] = omtauinv * (*f_pre)[i] + tauinv * f_eq;
    }
}

void
Node::stream(Lattice& lattice)
{
  // direction numbering scheme
  // 6 2 5
  // 3 0 1
  // 7 4 8

  unsigned int x = coord[0];
  unsigned int y = coord[1];
  unsigned int x_forward = 0;
  unsigned int y_forward = 0;

  (*f_adj)[0] = (*f_post)[0];
  
  for(unsigned int i = 1; i<dir; ++i){
    if(!bounce_back_dir[i]){
      x_forward = x + coeff[i][0];
      y_forward = y + coeff[i][1];
      lattice.get_node(x_forward, y_forward).set_f_adj(i, (*f_post)[i]);
    }
  }
}

void
Node::apply_BCs(const Lattice &lattice)
{
  unsigned int x_forward = 0;
  unsigned int y_forward = 0;
  for(unsigned int i = 0; i<dir; ++i){
    if(bounce_back_dir[i]){
      x_forward = coord[0] + coeff[i][0];
      y_forward = coord[1] + coeff[i][1];

      NodeType type = lattice.get_node(x_forward, y_forward).get_node_type();
      
      if(type == NodeType::solid ||
        type == NodeType::obstacle ||
        type == NodeType::inlet){
        apply_IBB(lattice, i);
        // apply_BB(lattice, i);
      }
      else if(type == NodeType::outlet){
        apply_ZouHe(lattice, i);// TODO: zouhe has to be outside the for loop no need for is
      }
      else{
        throw std::runtime_error("Invalid BCs type");
      }
    }
  }
}

void
Node::apply_IBB(const Lattice &lattice, unsigned int i)
{
  
  // Interpolated Bounce-Back 
  // take the velocity at the wall node
  unsigned int x_forward = coord[0] + coeff[i][0];
  unsigned int y_forward = coord[1] + coeff[i][1];
  double ux_wall = lattice.get_node(x_forward, y_forward).get_ux();
  double uy_wall = lattice.get_node(x_forward, y_forward).get_uy();

  NodeType type = lattice.get_node(x_forward, y_forward).get_node_type();
  if(type == NodeType::outlet){
    unsigned int x_backward = coord[0] + coeff[bb_indexes[i]][0];
    unsigned int y_backward = coord[1] + coeff[bb_indexes[i]][1];
    double ux_fluid = lattice.get_node(x_backward, y_backward).get_ux();
    double uy_fluid = lattice.get_node(x_backward, y_backward).get_uy();

    ux_wall = 0.5 * ux;// - 0.5 * ux_fluid;
    uy_wall = 0.5 * uy;// - 0.5 * uy_fluid; 
  }

  if(check_backward(lattice, coord[0], coord[1], i))
  {
    // since we have already collided and streamed, we take the post-collision value from f_adj
    // f_adj_post_coll = lattice.get_node(x_forward, y_forward).get_f_post(i);
    double f_adj_post_coll = (*f_adj)[i]; 
    
    (*f_adj)[bb_indexes[i]] = (2 * bounce_back_delta[i] * (*f_post)[i] + 
                    (1 - 2 * bounce_back_delta[i]) * f_adj_post_coll) * 
                    (bounce_back_delta[i] < 0.5) +
                    (1. / (2 * bounce_back_delta[i]) * (*f_post)[i] + 
                    ((2 * bounce_back_delta[i] - 1.) / (2 * bounce_back_delta[i])) * (*f_post)[bb_indexes[i]]) *
                    (bounce_back_delta[i] >= 0.5) - 
                    (ux_wall * coeff[i][0] + uy_wall * coeff[i][1]) * weights[i] * 6; // Wall velocity term (rho)
  }
  else
  {
    (*f_adj)[bb_indexes[i]] = (*f_post)[i];
  }
}

void
Node::apply_BB(const Lattice &lattice, unsigned int i)
{
  // Simple Bounce-Back
  unsigned int x_forward = coord[0] + coeff[i][0];
  unsigned int y_forward = coord[1] + coeff[i][1];
  double ux_wall = lattice.get_node(x_forward, y_forward).get_ux();
  double uy_wall = lattice.get_node(x_forward, y_forward).get_uy();
  double rho_w = rho;

  NodeType type = lattice.get_node(x_forward, y_forward).get_node_type();
  if(type == NodeType::outlet){
    unsigned int x_backward = coord[0] + coeff[bb_indexes[i]][0];
    unsigned int y_backward = coord[1] + coeff[bb_indexes[i]][1];
    double ux_fluid = lattice.get_node(x_backward, y_backward).get_ux();
    double uy_fluid = lattice.get_node(x_backward, y_backward).get_uy();

    ux_wall = (ux + ux_fluid)/2.;
    uy_wall = (ux + uy_fluid)/2.;
    rho_w =  (2 * ((*f_post)[1] + (*f_post)[5] + (*f_post)[8]) + (*f_post)[0] + (*f_post)[2] + (*f_post)[4])/ (1 - ux_wall);
  }
  
  if(check_backward(lattice, coord[0], coord[1], i))
  {
    (*f_adj)[bb_indexes[i]] = (*f_post)[i] - 
                  (ux_wall * coeff[i][0] + uy_wall * coeff[i][1]) * rho * weights[i] * 6;
  }
  else
  {
    (*f_adj)[bb_indexes[i]] = (*f_post)[i];
  }
}

void 
Node::apply_anti_BB(const Lattice &lattice, unsigned int i)
{
  // Anti Bounce-Back for outlet nodes
  unsigned int x_forward = coord[0] + coeff[i][0];
  unsigned int y_forward = coord[1] + coeff[i][1];
  double ux_wall = lattice.get_node(x_forward, y_forward).get_ux();
  double uy_wall = lattice.get_node(x_forward, y_forward).get_uy();
  double rho_w = rho;

  NodeType type = lattice.get_node(x_forward, y_forward).get_node_type();
  if(type == NodeType::outlet){
    unsigned int x_backward = coord[0] + coeff[bb_indexes[i]][0];
    unsigned int y_backward = coord[1] + coeff[bb_indexes[i]][1];
    double ux_fluid = lattice.get_node(x_backward, y_backward).get_ux();
    double uy_fluid = lattice.get_node(x_backward, y_backward).get_uy();
    double rho_fluid = lattice.get_node(x_backward, y_backward).get_rho();
    double ux_fluid2 = lattice.get_node(x_backward + coeff[bb_indexes[i]][0], y_backward + coeff[bb_indexes[i]][1]).get_ux();
    double uy_fluid2 = lattice.get_node(x_backward + coeff[bb_indexes[i]][0], y_backward + coeff[bb_indexes[i]][1]).get_uy();


    ux_wall = 1.5 * ux - 0.5 * ux_fluid;
    uy_wall = 1.5 * uy - 0.5 * uy_fluid;
    rho_w = 1.5 * rho - 0.5 * rho_fluid;
  }
  
  if(check_backward(lattice, coord[0], coord[1], i))
  {
    (*f_adj)[bb_indexes[i]] = -(*f_post)[i]  +
                          2 * weights[i] * rho_w *
                          (1 + 4.5 * (coeff[i][0] * ux_wall + coeff[i][1] * uy_wall) * (coeff[i][0] * ux_wall + coeff[i][1] * uy_wall) -
                          3.5 * (ux_wall * ux_wall + uy_wall * uy_wall)); // The output of ux_out is fine ux is positive exiting the domain
  }
  else
  {
    (*f_adj)[bb_indexes[i]] = (*f_post)[i];
  }
}

//TODO: scrivere una funzione che determini per i nodi boundary la posizione dell'outlet:
//        - se sono ai corner o ai lati, fai switch e applica il caso specifico di zou-he 
void 
Node::apply_ZouHe(const Lattice &lattice, unsigned int i)
{
  double u_wall = 0.0;
  unsigned int x_backward = 0;
  unsigned int y_backward = 0;

  switch (zou_he_type)
  {
    case ZouHeType::none:
      // No Zou-He boundary condition
      break;
    case ZouHeType::right:
      rho = 1;
      u_wall = (*f_adj)[0] + (*f_adj)[2] + (*f_adj)[4] + 2.0 * ((*f_adj)[1] + (*f_adj)[5] + (*f_adj)[8]) - 1.0;
      (*f_adj)[3] = (*f_adj)[1] - 2.0 / 3.0 * ux_wall;
      (*f_adj)[6] = (*f_adj)[8] - 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * ux_wall;
      (*f_adj)[7] = (*f_adj)[5] + 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * ux_wall;
      break;
    case ZouHeType::top:
      rho = 1;
      u_wall = (*f_adj)[0] + (*f_adj)[1] + (*f_adj)[3] + 2.0 * ((*f_adj)[2] + (*f_adj)[5] + (*f_adj)[6]) - 1.0;
      (*f_adj)[4] = (*f_adj)[2] - 2.0 / 3.0 * u_wall;
      (*f_adj)[7] = (*f_adj)[5] - 0.5 * ((*f_adj)[3] - (*f_adj)[1]) - 1.0 / 6.0 * ux_wall;
      (*f_adj)[8] = (*f_adj)[6] + 0.5 * ((*f_adj)[3] - (*f_adj)[1]) - 1.0 / 6.0 * ux_wall;
      break;
      break;
    case ZouHeType::left:
      rho = 1;
      u_wall = (*f_adj)[0] + (*f_adj)[2] + (*f_adj)[4] + 2.0 * ((*f_adj)[3] + (*f_adj)[6] + (*f_adj)[7]) - 1.0;
      (*f_adj)[1] = (*f_adj)[3] - 2.0 / 3.0 * u_wall;
      (*f_adj)[5] = (*f_adj)[7] - 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * ux_wall;
      (*f_adj)[8] = (*f_adj)[6] + 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * ux_wall;
      break;
      break;
    case ZouHeType::bottom:
      rho = 1;
      u_wall = (*f_adj)[0] + (*f_adj)[1] + (*f_adj)[3] + 2.0 * ((*f_adj)[4] + (*f_adj)[7] + (*f_adj)[8]) - 1.0;
      (*f_adj)[2] = (*f_adj)[4] - 2.0 / 3.0 * u_wall;
      (*f_adj)[5] = (*f_adj)[7] - 0.5 * ((*f_adj)[1] - (*f_adj)[3]) - 1.0 / 6.0 * ux_wall;
      (*f_adj)[6] = (*f_adj)[8] + 0.5 * ((*f_adj)[1] - (*f_adj)[3]) - 1.0 / 6.0 * ux_wall;
      break;
    case ZouHeType::top_right:
        x_backward = coord[0] - 1;
        y_backward = coord[1] + 1;
        rho = 1;

        f.at() = f.at() - 2.0 / 3.0 * rho * macroU.at(0);
        f.at() = f.at() - 2.0 / 3.0 * rho * macroU.at(1);
        f.at() = f.at() - 1.0 / 6.0 * rho * macroU.at(0) - 1.0 / 6.0 * rho * macroU.at(1);
        f.at() = 0;
        f.at() = 0;
        f.at() = rho - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
      break;
    case ZouHeType::top_left:
      // Apply Zou-He boundary condition for top-left corner
      break;
    case ZouHeType::bottom_left:
      // Apply Zou-He boundary condition for bottom-left corner
      break;
    case ZouHeType::bottom_right:
      // Apply Zou-He boundary condition for bottom-right corner
      break;
    default:
      throw std::runtime_error("Invalid ZouHeType");
  }



  if(check_backward(lattice, coord[0], coord[1], i))
  {

    rho = 1;
    ux_wall = (*f_adj)[0] + (*f_adj)[2] + (*f_adj)[4] + 2.0 * ((*f_adj)[1] + (*f_adj)[5] + (*f_adj)[8]) - 1.0;
    (*f_adj)[3] = (*f_adj)[1] - 2.0 / 3.0 * ux_wall;
    (*f_adj)[6] = (*f_adj)[8] - 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * ux_wall;
    (*f_adj)[7] = (*f_adj)[5] + 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * ux_wall;
  }
  else
  {
    (*f_adj)[bb_indexes[i]] = (*f_post)[i];
  }
}

//TODO: check if necessary
bool
Node::check_backward(const Lattice &lattice, unsigned int x,unsigned int y, unsigned int i)
{
  unsigned int x_backward = coord[0] + coeff[bb_indexes[i]][0];
  unsigned int y_backward = coord[1] + coeff[bb_indexes[i]][1];

  NodeType type = lattice.get_node(x_backward, y_backward).get_node_type();

  if(type == NodeType::fluid || type == NodeType::boundary)
    return true;
  else
    return false;
}

void
Node::update_f()
{
  f_pre.swap(f_adj);
}

void
Node::compute_physical_quantities()
{
  // compute rho and U usign the equilibrium distribution
  double rho_ = 0.0;
  for (unsigned int i = 0; i < dir; ++i)
  {
      rho_ += (*f_pre)[i];
  }
  rho = rho_;
  double rhoinv = 1.0/rho;

  double u = 0.;
  double v = 0.;
  for (unsigned int i = 0; i < dir; ++i)
  {
      u += rhoinv*coeff[i][0]*(*f_pre)[i];
      v += rhoinv*coeff[i][1]*(*f_pre)[i];
  }

  ux = u;
  uy = v;
}

void
Node::compute_drag_and_lift(const Lattice &lattice)
{
  // compute drag and lift
  drag = 0.0;
  lift = 0.0;

  unsigned int x_forward = 0;
  unsigned int y_forward = 0;

  for(unsigned int i = 0; i<dir; ++i){
    x_forward = coord[0] + coeff[i][0];
    y_forward = coord[1] + coeff[i][1];
    if(lattice.get_node(x_forward, y_forward).get_node_type() == NodeType::obstacle){
      drag += coeff[i][0] * (*f_pre)[i] - coeff[bb_indexes[i]][0] * (*f_adj)[bb_indexes[i]];
      lift -= coeff[i][1] * (*f_pre)[i] - coeff[bb_indexes[i]][1] * (*f_adj)[bb_indexes[i]];
    }
  }
}

void
Node::set_bounce_back_properties(std::vector<bool> bounce_back_dir_, 
                                  std::vector<double> bounce_back_delta_)
{
  bounce_back_dir = bounce_back_dir_;
  bounce_back_delta = bounce_back_delta_;
  set_zou_he_type();
}

void
Node::set_zou_he_type()
{
  if(bounce_back_dir[5] && bounce_back_dir[1] && bounce_back_dir[8])
    zou_he_type = ZouHeType::right;
  else if(bounce_back_dir[6] && bounce_back_dir[2] && bounce_back_dir[5])
    zou_he_type = ZouHeType::top;
  else if(bounce_back_dir[7] && bounce_back_dir[3] && bounce_back_dir[6])
    zou_he_type = ZouHeType::left;
  else if(bounce_back_dir[8] && bounce_back_dir[4] && bounce_back_dir[7])
    zou_he_type = ZouHeType::bottom;
  else if(bounce_back_dir[2] && bounce_back_dir[5] && bounce_back_dir[1])
    zou_he_type = ZouHeType::top_right;
  else if(bounce_back_dir[2] && bounce_back_dir[6] && bounce_back_dir[3])
    zou_he_type = ZouHeType::top_left;
  else if(bounce_back_dir[4] && bounce_back_dir[7] && bounce_back_dir[3])
    zou_he_type = ZouHeType::bottom_left;
  else if(bounce_back_dir[4] && bounce_back_dir[8] && bounce_back_dir[1])
    zou_he_type = ZouHeType::bottom_right;
  else
    zou_he_type = ZouHeType::none;
}