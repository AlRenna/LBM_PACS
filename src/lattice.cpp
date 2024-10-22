/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "lattice.hpp"

Lattice::Lattice(unsigned int nx_, unsigned int ny_,
                 double nu_)
  : nx(nx_),
    ny(ny_),
    nu(nu_),
    tau(3.0*nu+0.5)
{
  nodes.resize(nx*ny);
  node_types.resize(nx*ny, NodeType::solid);
  boundary_node_dir.resize(nx*ny, std::vector<bool>(9, false));
  ux_in.resize(nx*ny, 0.);
  uy_in.resize(nx*ny, 0.);
  rho_in.resize(nx*ny, 1.0);
  ux_out.resize(nx*ny, 0.);
  uy_out.resize(nx*ny, 0.);
  rho_out.resize(nx*ny, 1.0);
}

Lattice::set_ICs_&_BCs(std::vector<double> ux_in_, 
                            std::vector<double> uy_in_, 
                            std::vector<double> rho_in_,
                            std::vector<NodeType> node_types_,
                            std::vector<std::vector<bool>> boundary_node_dir_)
{
  ux_in = ux_in_;
  uy_in = uy_in_;
  rho_in = rho_in_;
  node_types = node_types_;
  boundary_node_dir = boundary_node_dir_;
}

Lattice::populate_Nodes()
{
  for(unsigned int y = 0; y<ny; ++y){
    for(unsigned int x = 0; x<nx; ++x){
      unsigned int index = scalar_index(x, y);
      nodes[index] = Node(node_types[index], {x, y}, ux[index], uy[index], rho[index]);
      nodes[index].set_boundary_node_dir(boundary_node_dir[index]);
      nodes[index].init();
    }
  }
}

Lattice::run()
{
  for(/*time loop*/)
  {
    for(unsigned int y = 0; y<ny; ++y)
    {
      for(unsigned int x = 0; x<nx; ++x)
      {
        unsigned int index = scalar_index(x, y);
        if(node_types[index] != NodeType::solid)
        {
          
          // TODO: check order of operations
          nodes[index].apply_bc(); // inlet e BCs (zou he)
          nodes[index].compute_physical_quantities();
          nodes[index].load_adjacent_velocity_distributions(*this);
          nodes[index].collide_stream(*this);

          // nodes[index].compute_integrals();
          // nodes[index].save(*this);

          
        }
      }
    }
  }
}