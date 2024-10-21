#include "node.hpp"

Node::Node(NodeType node_type_, std::vector<int> coord_,
        double ux_, double uy_, double rho_)
  : node_type(node_type_),
    coord(coord_),
    ux(ux_),
    uy(uy_),
    rho(rho_)
{
  f.resize(dir, 0.);
  f_adj.resize(dir, 0.);
}

Node::init()
{
  // initialize f

  std::vector<double> w_rho = weights*rho;
  double omusq = 1.0 - 1.5*(ux*ux+uy*uy);
  
  double tu = 3.0*ux;
  double tv = 3.0*uy;

  for(unsigned int i = 0; i<dir; ++i){
    f[i] = w_rho[i]*(omusq + 
            (coeff[i][0]*tu + coeff[i][1]*tv)*(1.0+0.5*(coeff[i][0]*tu + coeff[i][1]*tv)));
  }
}

Node::load_adjacent_velocity_distributions(const Lattice &lattice)
{
  nx = lattice.get_nx();
  ny = lattice.get_ny();
  x = coord[0];
  y = coord[1];

  unsigned int x_dx = (x+1)%nx;
  unsigned int y_down = (y+1)%ny;
  unsigned int x_sx = (nx+x-1)%nx;
  unsigned int y_up = (ny+y-1)%ny;
  
  // direction numbering scheme
  // 6 2 5
  // 3 0 1
  // 7 4 8
  
  f_adj[0] = f[0];

  // load populations from adjacent nodes  
  f_adj[1] = lattice.get_node(x_sx, y     ).get_f(1);
  f_adj[2] = lattice.get_node(x,    y_up  ).get_f(2);
  f_adj[3] = lattice.get_node(x_dx, y     ).get_f(3);
  f_adj[4] = lattice.get_node(x,    y_down).get_f(4);
  f_adj[5] = lattice.get_node(x_sx, y_up  ).get_f(5);
  f_adj[6] = lattice.get_node(x_dx, y_up  ).get_f(6);
  f_adj[7] = lattice.get_node(x_dx, y_down).get_f(7);
  f_adj[8] = lattice.get_node(x_sx, y_down).get_f(8);
}

Node::compute_physical_quantities()
{
  // compute rho and U usign the equilibrium distribution
  double rho_ = 0.0;
  for (unsigned int i = 0; i < 9; ++i)
  {
      rho_ += f_adj[i];
  }
  rho = rho_;
  double rhoinv = 1.0/rho;

  double u = 0.;
  double v = 0.;
  for (unsigned int i = 0; i < 9; ++i)
  {
      u += rhoinv*coeff[i][0]*f_adj[i];
      v += rhoinv*coeff[i][1]*f_adj[i];
  }

  ux = u;
  uy = v;
}

Node::collide_stream(const Lattice &lattice)
{
  double tau = lattice.get_tau();
  const double tauinv = 1./tau; // 1/tau
  const double omtauinv = 1.0-tauinv;     // 1 - 1/tau
  
  // Collision and streaming step: now compute and relax to equilibrium 
  // BGK collision operator
  // f_i(t+dt) = f_i(t) - 1/tau (f_i(t) - f_i^eq(t))

  // temporary variables
  std::vector<double> w_tau_rho = weights*rho_*tauinv;
  double omusq = 1.0 - 1.5*(u*u+v*v); // 1-(3/2)u.u
  
  double tu = 3.0*u;
  double tv = 3.0*v;
    
  for(unsigned int i = 0; i<9; ++i){
    f[i] = omtauinv*f_adj[i] +
            w_tau_rho[i]*(omusq + 
            (coeff[i][0]*tu + coeff[i][1]*tv)*(1.0+0.5*(coeff[i][0]*tu + coeff[i][1]*tv)));
  }
}