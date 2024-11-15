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
    rho(rho_)
{
  f = std::make_unique<std::vector<double>>(dir, 0.);
  f_adj = std::make_unique<std::vector<double>>(dir, 0.);
  boundary_node_dir.resize(dir, false);
  boundary_node_delta.resize(dir, 0.);
}

void
Node::init_equilibrium()
{
  // initialize f to equilibrium distribution (only at step 0)
  // f_i = w_i * rho * (1  + 3 * (c_i.u) + 4.5 * (c_i.u)^2 - 1.5 * u^2)

  for(unsigned int i = 0; i<dir; ++i){
    (*f)[i] = weights[i] * rho * (1. + 3.0 * (coeff[i][0] * ux + coeff[i][1] * uy) + 
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
  for(unsigned int i = 0; i<9; ++i){
    f_eq =  weights[i] * rho * (1. + 3.0 * (coeff[i][0] * ux + coeff[i][1] * uy) + 
            4.5 * (coeff[i][0] * ux + coeff[i][1] * uy) * (coeff[i][0] * ux + coeff[i][1] * uy) - 
            1.5 * (ux * ux + uy * uy));
    (*f)[i] = omtauinv * (*f)[i] + tauinv * f_eq;
    }
}

void
Node::apply_IBB(const Lattice &lattice)
{
  double f_adj_post_coll = 0.0;
  unsigned int x_adg = 0;
  unsigned int y_adg = 0;
  double ux_wall = 0.;
  double uy_wall = 0.;

  // Interpolated Bounce-Back
  for(unsigned int i = 1; i<dir; ++i){
    if(boundary_node_dir[i]){

      // take the velocity at the wall node
      x_adg = coord[0] + coeff[i][0];
      y_adg = coord[1] + coeff[i][1];
      ux_wall = lattice.get_node(x_adg, y_adg).get_ux();
      uy_wall = lattice.get_node(x_adg, y_adg).get_uy();

      // since we have already collided and streamed, we take the post-collision value from f_adj
      // f_adj_post_coll = lattice.get_node(x_adg, y_adg).get_f(i);
      f_adj_post_coll = (*f_adj)[i]; 
      
      (*f_adj)[bb_indexes[i]] = (2 * boundary_node_delta[i] * (*f)[i] + 
                      (1 - 2 * boundary_node_delta[i]) * f_adj_post_coll) * 
                      (boundary_node_delta[i] < 0.5) +
                      (1. / (2 * boundary_node_delta[i]) * (*f)[i] + 
                      ((2 * boundary_node_delta[i] - 1.) / (2 * boundary_node_delta[i])) * (*f)[bb_indexes[i]]) *
                      (boundary_node_delta[i] >= 0.5) - 
                      (ux_wall * coeff[i][0] + uy_wall * coeff[i][1]) * weights[i] * 6; // Wall velocity term (rho)
    }
  }
}

void
Node::apply_BB(const Lattice &lattice)
{
  unsigned int x_adg = 0;
  unsigned int y_adg = 0;
  double ux_wall = 0.;
  double uy_wall = 0.;

  for(unsigned int i = 1; i<dir; ++i){
    if(boundary_node_dir[i]){
      x_adg = coord[0] + coeff[i][0];
      y_adg = coord[1] + coeff[i][1];
      ux_wall = lattice.get_node(x_adg, y_adg).get_ux();
      uy_wall = lattice.get_node(x_adg, y_adg).get_uy();

      (*f_adj)[bb_indexes[i]] = (*f)[i] - 
                      (ux_wall * coeff[i][0] + uy_wall * coeff[i][1]) * rho * weights[i] * 6;
    }
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
  unsigned int x_adg = 0;
  unsigned int y_adg = 0;

  (*f_adj)[0] = (*f)[0];

  for(unsigned int i = 1; i<dir; ++i){
    if(!boundary_node_dir[i]){
      x_adg = x + coeff[i][0];
      y_adg = y + coeff[i][1];
      lattice.get_node(x_adg, y_adg).set_f_adj(i, (*f)[i]);
      // std::cout << "Streaming: " << x << " " << y << " in " << x_adg << " " << y_adg << " :" << (*f)[i] << std::endl;
    }
  }
}

void
Node::compute_physical_quantities()
{
  // compute rho and U usign the equilibrium distribution
  double rho_ = 0.0;
  for (unsigned int i = 0; i < 9; ++i)
  {
      rho_ += (*f)[i];
  }
  rho = rho_;
  double rhoinv = 1.0/rho;

  double u = 0.;
  double v = 0.;
  for (unsigned int i = 0; i < 9; ++i)
  {
      u += rhoinv*coeff[i][0]*(*f)[i];
      v += rhoinv*coeff[i][1]*(*f)[i];
  }

  ux = u;
  uy = v;
}

void
Node::update_f()
{
  f.swap(f_adj);
}

void
Node::set_boundary_node_properties(std::vector<bool> boundary_node_dir_, 
                                  std::vector<double> boundary_node_delta_,
                                  double dt)
{
  boundary_node_dir = boundary_node_dir_;
  boundary_node_delta = boundary_node_delta_;
}
