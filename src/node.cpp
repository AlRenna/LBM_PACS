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
    zou_he_type(ZouHeType::none)
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
        continue; // ZouHe is applied after the loop for outlet nodes
      }
      else{
        throw std::runtime_error("Invalid BCs type");
      }
    }
  }

  apply_ZouHe(lattice);
}

void
Node::apply_IBB(const Lattice &lattice, unsigned int i)
{
  
  // Interpolated Bounce-Back 
  // take the velocity at the wall node
  unsigned int x_forward = coord[0] + coeff[i][0];
  unsigned int y_forward = coord[1] + coeff[i][1];
  double current_time = lattice.get_time();
  double ux_wall = lattice.get_node(x_forward, y_forward).get_ux() * (1/ (1 + std::exp(-25 *(current_time - 0.2))));
  double uy_wall = lattice.get_node(x_forward, y_forward).get_uy() * (1/ (1 + std::exp(-25 *(current_time - 0.2))));

  // check if the node in the backward direction is a fluid or boundary node
  if(!bounce_back_dir[bb_indexes[i]])
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
  double current_time = lattice.get_time();
  double ux_wall = lattice.get_node(x_forward, y_forward).get_ux() * (1/ (1 + std::exp(-25 *(current_time - 0.2))));
  double uy_wall = lattice.get_node(x_forward, y_forward).get_uy() * (1/ (1 + std::exp(-25 *(current_time - 0.2))));
  
  // check if the node in the backward direction is a fluid or boundary node
  if(!bounce_back_dir[bb_indexes[i]])
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
Node::apply_ZouHe(const Lattice &lattice)
{
  double u_wall = 0.0;
  int x = coord[0];
  int y = coord[1];
  double rho_wall = 0.0;

  switch (zou_he_type)
  {
    case ZouHeType::none:
      // No Zou-He boundary condition
      break;
    case ZouHeType::right:
      rho_wall = lattice.get_node(x + coeff[1][0], y + coeff[1][1]).get_rho();
      u_wall = (*f_adj)[0] + (*f_adj)[2] + (*f_adj)[4] + 2.0 * ((*f_adj)[1] + (*f_adj)[5] + (*f_adj)[8]) - rho_wall;
      (*f_adj)[3] = (*f_adj)[1] - 2.0 / 3.0 * u_wall;
      (*f_adj)[6] = (*f_adj)[8] - 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * u_wall;
      (*f_adj)[7] = (*f_adj)[5] + 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * u_wall;
      break;
    case ZouHeType::top:
      rho_wall = lattice.get_node(x + coeff[2][0], y + coeff[2][1]).get_rho();
      u_wall = (*f_adj)[0] + (*f_adj)[1] + (*f_adj)[3] + 2.0 * ((*f_adj)[2] + (*f_adj)[5] + (*f_adj)[6]) - rho_wall;
      (*f_adj)[4] = (*f_adj)[2] - 2.0 / 3.0 * u_wall;
      (*f_adj)[7] = (*f_adj)[5] - 0.5 * ((*f_adj)[3] - (*f_adj)[1]) - 1.0 / 6.0 * u_wall;
      (*f_adj)[8] = (*f_adj)[6] + 0.5 * ((*f_adj)[3] - (*f_adj)[1]) - 1.0 / 6.0 * u_wall;
      break;
    case ZouHeType::left:
      rho_wall = lattice.get_node(x + coeff[3][0], y + coeff[3][1]).get_rho();
      u_wall = (*f_adj)[0] + (*f_adj)[2] + (*f_adj)[4] + 2.0 * ((*f_adj)[3] + (*f_adj)[6] + (*f_adj)[7]) - rho_wall;
      (*f_adj)[1] = (*f_adj)[3] - 2.0 / 3.0 * u_wall;
      (*f_adj)[5] = (*f_adj)[7] - 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * u_wall;
      (*f_adj)[8] = (*f_adj)[6] + 0.5 * ((*f_adj)[2] - (*f_adj)[4]) - 1.0 / 6.0 * u_wall;
      break;
    case ZouHeType::bottom:
      rho_wall = lattice.get_node(x + coeff[4][0], y + coeff[4][1]).get_rho();
      u_wall = (*f_adj)[0] + (*f_adj)[1] + (*f_adj)[3] + 2.0 * ((*f_adj)[4] + (*f_adj)[7] + (*f_adj)[8]) - rho_wall;
      (*f_adj)[2] = (*f_adj)[4] - 2.0 / 3.0 * u_wall;
      (*f_adj)[5] = (*f_adj)[7] - 0.5 * ((*f_adj)[1] - (*f_adj)[3]) - 1.0 / 6.0 * u_wall;
      (*f_adj)[6] = (*f_adj)[8] + 0.5 * ((*f_adj)[1] - (*f_adj)[3]) - 1.0 / 6.0 * u_wall;
      break;
    default:
      throw std::runtime_error("Invalid ZouHeType");
  }

}

void
Node::update_f()
{
  f_pre.swap(f_adj);
}

void
Node::compute_physical_quantities()
{
  // compute rho and U using the equilibrium distribution
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
}

void
Node::set_zou_he_type(Lattice &lattice)
{
  // if(bounce_back_dir[5] && bounce_back_dir[1] && bounce_back_dir[8])
  //   zou_he_type = ZouHeType::right;
  // else if(bounce_back_dir[6] && bounce_back_dir[2] && bounce_back_dir[5])
  //   zou_he_type = ZouHeType::top;
  // else if(bounce_back_dir[7] && bounce_back_dir[3] && bounce_back_dir[6])
  //   zou_he_type = ZouHeType::left;
  // else if(bounce_back_dir[8] && bounce_back_dir[4] && bounce_back_dir[7])
  //   zou_he_type = ZouHeType::bottom;
  // else if(bounce_back_dir[2] && bounce_back_dir[5] && bounce_back_dir[1])
  //   zou_he_type = ZouHeType::top_right;
  // else if(bounce_back_dir[2] && bounce_back_dir[6] && bounce_back_dir[3])
  //   zou_he_type = ZouHeType::top_left;
  // else if(bounce_back_dir[4] && bounce_back_dir[7] && bounce_back_dir[3])
  //   zou_he_type = ZouHeType::bottom_left;
  // else if(bounce_back_dir[4] && bounce_back_dir[8] && bounce_back_dir[1])
  //   zou_he_type = ZouHeType::bottom_right;
  // else
  //   zou_he_type = ZouHeType::none;
  int x = coord[0];
  int y = coord[1];
  if(bounce_back_dir[1] && lattice.get_node(x + coeff[1][0], y + coeff[1][1]).get_node_type() == NodeType::outlet)
    zou_he_type = ZouHeType::right;
  else if(bounce_back_dir[2] && lattice.get_node(x + coeff[2][0], y + coeff[2][1]).get_node_type() == NodeType::outlet)
    zou_he_type = ZouHeType::top;
  else if(bounce_back_dir[3] && lattice.get_node(x + coeff[3][0], y + coeff[3][1]).get_node_type() == NodeType::outlet)   
    zou_he_type = ZouHeType::left;
  else if(bounce_back_dir[4] && lattice.get_node(x + coeff[4][0], y + coeff[4][1]).get_node_type() == NodeType::outlet)
    zou_he_type = ZouHeType::bottom;
  else
    zou_he_type = ZouHeType::none;
}