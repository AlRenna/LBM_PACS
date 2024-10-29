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
  boundary_node_dir.resize(nx*ny, std::vector<bool>(Node::dir, false));
  ux_in.resize(nx*ny, 0.);
  uy_in.resize(nx*ny, 0.);
  rho_in.resize(nx*ny, 1.0);
  ux_out.resize(nx*ny, 0.);
  uy_out.resize(nx*ny, 0.);
  rho_out.resize(nx*ny, 1.0);
}

void
Lattice::load_ICs_and_BCs(const std::vector<double>& ux_in_, 
                          const std::vector<double>& uy_in_, 
                          const std::vector<double>& rho_in_,
                          const std::string& filename_nodes)
{
  ux_in = ux_in_;
  uy_in = uy_in_;
  rho_in = rho_in_;
  readNodesFromCSV(filename_nodes);

  // Populate the nodes
  populate_Nodes();
}

void 
Node::readNodesFromCSV(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
      throw std::runtime_error("Could not open file");
  }

  std::string line;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string token;

    // Read coordinates
    std::getline(ss, token, ',');
    unsigned int x = std::stoi(token);
    std::getline(ss, token, ',');
    unsigned int y = std::stoi(token);

    // Read NodeType
    std::getline(ss, token, ',');
    NodeType nodeType = static_cast<NodeType>(std::stoi(token));

    // Read the last 8 columns
    for (int i = 1; i < Node::dir; ++i) {
      std::getline(ss, token, ',');
      double value = std::stod(token);
      boundary_node_delta[i] = value;
      boundary_node_dir[i] = (value > 0.0);
    }
  }

  file.close();
}

void
Lattice::populate_Nodes()
{
  for(unsigned int y = 0; y<ny; ++y){
    for(unsigned int x = 0; x<nx; ++x){
      unsigned int index = scalar_index(x, y);
      nodes[index] = Node(node_types[index], {x, y}, ux_in[index], uy_in[index], rho_in[index]);
      nodes[index].set_boundary_node_properties(boundary_node_dir[index], boundary_node_delta[index], dt);
      nodes[index].init_equilibrium();
    }
  }
}

void
Lattice::run()
{
  unsigned int iter = 0;
  while(iter < max_iter)
  {
    for(unsigned int y = 0; y<ny; ++y)
    {
      for(unsigned int x = 0; x<nx; ++x)
      {
        unsigned int index = scalar_index(x, y);
        if(node_types[index] != NodeType::solid)
        {
          
          // nodes[index].apply_bc(); // inlet e BCs (zou he)
          nodes[index].collide(*this);
          nodes[index].stream(*this);
        }
      }
    }

    // We separate the streaming step from the collision step
    // We need the updated information on all the node to procede
    for(unsigned int y = 0; y<ny; ++y)
    {
      for(unsigned int x = 0; x<nx; ++x)
      {
        unsigned int index = scalar_index(x, y);
        if(node_types[index] != NodeType::solid)
        {
          if(node_types[index] == NodeType::boundary)
            nodes[index].apply_IBB(); // BCs (zou he)
          
          nodes[index].compute_physical_quantities();
          // nodes[index].compute_integrals();

          // nodes[index].save(*this);
          nodes[index].update_f();          
        }
      }
    }

    iter = iter + 1;
  }
}
