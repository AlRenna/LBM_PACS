/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "lattice.hpp"

//TODO: cambiare i params.json in modo che contenga anche il dt e il max_iter, salvataggio
//TODO: vedere inlet outlet
//TODO: Fare un esempio con ostacolo
// TODO: calcolo lift e drag
Lattice::Lattice(unsigned int nx_, unsigned int ny_,
                 double nu_)
  : nx(nx_),
    ny(ny_),
    nu(nu_),
    tau(3.0*nu+0.5)
{
  nodes.resize(nx*ny);
  node_types.resize(nx*ny, NodeType::solid);
  boundary_node_delta.resize(nx*ny, std::vector<double>(Node::dir, 0.0));
  boundary_node_dir.resize(nx*ny, std::vector<bool>(Node::dir, false));
  ux_in.resize(nx*ny, 0.);
  uy_in.resize(nx*ny, 0.);
  rho_in.resize(nx*ny, 1.0);
  ux_out.resize(nx*ny, 0.);
  uy_out.resize(nx*ny, 0.);
  rho_out.resize(nx*ny, 1.0);
}

Lattice::Lattice()
{
  std::ifstream param_file("params.json");
  if (!param_file.is_open()) {
    throw std::runtime_error("Could not open params.json");
  }

  nlohmann::json param_json;
  param_file >> param_json;

  if (param_json.find("new_nx") == param_json.end() || 
      param_json.find("new_ny") == param_json.end() || 
      param_json.find("nu") == param_json.end()) {
    throw std::runtime_error("param.json does not contain required parameters");
  }

  nx = param_json["new_nx"];
  ny = param_json["new_ny"];
  nu = param_json["nu"];
  dt = param_json["dt"];
  max_iter = param_json["max_iter"];
  T_final = max_iter * dt;
  tau = 3.0 * nu + 0.5;
  nodes.resize(nx*ny);
  node_types.resize(nx*ny, NodeType::solid);
  boundary_node_delta.resize(nx*ny, std::vector<double>(Node::dir, 0.0));
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
  ux_out = ux_in;
  uy_out = uy_in;
  rho_out = rho_in;
  readNodesFromCSV(filename_nodes);

  // Populate the nodes
  populate_Nodes();
}

void 
Lattice::readNodesFromCSV(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
      throw std::runtime_error("Could not open file");
  }

  std::string line;
  // Skip the first line
  std::getline(file, line);

  while (std::getline(file, line)) {
    if (line.empty()) {
      break;
    }
    std::stringstream ss(line);
    std::string token;

    // Read coordinates
    std::getline(ss, token, ',');
    unsigned int x = std::stoi(token);
    std::getline(ss, token, ',');
    unsigned int y = std::stoi(token);
    unsigned int index = scalar_index(x, y);

    // Read NodeType
    std::getline(ss, token, ',');
    NodeType nodeType = static_cast<NodeType>(std::stoi(token));
    node_types[index] = nodeType;

    // Read the last 8 columns
    for (int i = 1; i < Node::dir; ++i) {
      std::getline(ss, token, ',');
      double value = std::stod(token);
      boundary_node_delta[index][i] = value;
      boundary_node_dir[index][i] = (value > std::numeric_limits<double>::epsilon());
    }
  }

  file.close();
}

void
Lattice::populate_Nodes()
{
  std::cout << "Populating Nodes" << std::endl;
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
  std::cout << "Running simulation\n" << std::endl;
  unsigned int iter = 0;

  // Delete the output_results directory if it exists
  if (std::filesystem::exists("output_results")) {
    std::filesystem::remove_all("output_results");
    std::filesystem::create_directory("output_results");
  }
  else{
    std::filesystem::create_directory("output_results");
  }
  
  std::string u_filename = "output_results/velocity_out.txt";
  std::string rho_filename = "output_results/rho_out.txt";

  std::ofstream u_file(u_filename);
  std::ofstream rho_file(rho_filename);

  // Save the initial conditions
  writeResults(u_file, rho_file);
  //writeResults(iter);

  iter = iter + 1;
  
  while(iter <= max_iter)
  {
    std::cout << "Iteration: " << iter << std::endl;
    std::cout << "Time: " << iter*dt << std::endl;
    std::cout << "Collision and streaming" << std::endl;
    for(unsigned int y = 0; y<ny; ++y)
    {
      for(unsigned int x = 0; x<nx; ++x)
      {
        unsigned int index = scalar_index(x, y);
        if(node_types[index] != NodeType::solid)
        {
          nodes[index].collide(*this);
          // if(node_types[index] == NodeType::boundary)
          // {
          //   nodes[index].apply_BB(*this);
          // }
          nodes[index].stream(*this);
        }
      }
    }

    std::cout << "Physical quantities evaluation\n" << std::endl;
    // We separate the streaming step from the collision step
    // We need the updated information on all the node to procede
    for(unsigned int y = 0; y<ny; ++y)
    {
      for(unsigned int x = 0; x<nx; ++x)
      {
        unsigned int index = scalar_index(x, y);
        if(node_types[index] != NodeType::solid)
        {
          if(node_types[index] == NodeType::boundary){
            nodes[index].apply_IBB(*this);
          }
          nodes[index].update_f();

          nodes[index].compute_physical_quantities();
          // nodes[index].compute_integrals();

          ux_out[index] = nodes[index].get_ux();
          uy_out[index] = nodes[index].get_uy();
          rho_out[index] = nodes[index].get_rho();

                    
        }
      }
    }

    // save the results every 5 iterations
    if( iter%20 == 0 || iter == max_iter-1)
    {
      std::cout << "Writing results" << std::endl;
      writeResults(u_file, rho_file);
    }
    iter = iter + 1;
  }
}

void 
Lattice::writeResults(std::ofstream &file_u, std::ofstream &file_rho) {
  // Save ux_out
  for (unsigned int y = 0; y < ny; ++y) {
    for (unsigned int x = 0; x < nx; ++x) {
      unsigned int index = scalar_index(x, y);
      file_u << std::sqrt(ux_out[index] * ux_out[index] + uy_out[index] * uy_out[index]) << " ";
      file_rho << rho_out[index] << " ";
    }
  }
  file_u << "\n";
  file_rho << "\n";
  
}
