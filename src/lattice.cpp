/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "lattice.hpp"

// TODO: fare un nuovo costruttore che prende filename per geometria lattice, condizioni iniziali e al bordo
// TODO: rivedere funzione per il salvataggio dei risultati
// TODO: scrivi script per la lettura dei file .csv da python e crea un animazione

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
      boundary_node_dir[index][i] = (value > 0.0);
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
      // TODO: capire come inizializzare il file txt delle condizioni initiali considerando anche i nodi solidi
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
  while(iter < max_iter)
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
          
          // nodes[index].apply_bc(); // inlet e BCs (zou he)
          nodes[index].collide(*this);
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
          if(node_types[index] == NodeType::boundary)
            nodes[index].apply_IBB(*this); // BCs (zou he)
          
          nodes[index].compute_physical_quantities();
          // nodes[index].compute_integrals();

          // nodes[index].save(*this);
          // TODO: scrivere nelle variabili di output le quantità fisiche
          ux_out[index] = nodes[index].get_ux();
          uy_out[index] = nodes[index].get_uy();
          rho_out[index] = nodes[index].get_rho();

          nodes[index].update_f();          
        }
      }
    }

    std::cout << "Writing results" << std::endl;
    if( iter%5 == 0 || iter == max_iter-1)
    {
    // save the 
      writeResults(iter);
    }
    iter = iter + 1;
  }
}

void 
Lattice::writeResults(const unsigned int iter) {
  // Create directory if it doesn't exist
  std::filesystem::create_directory("output_results");

  double curr_time_step = iter * dt;

  // Create filenames based on the iteration
  std::string ux_filename = "output_results/ux_out_" + std::to_string(curr_time_step) + ".csv";
  std::string uy_filename = "output_results/uy_out_" + std::to_string(curr_time_step) + ".csv";
  std::string rho_filename = "output_results/rho_out_" + std::to_string(curr_time_step) + ".csv";

  // Save ux_out
  std::ofstream ux_file(ux_filename);
  for (unsigned int y = 0; y < ny; ++y) {
    for (unsigned int x = 0; x < nx; ++x) {
      ux_file << ux_out[scalar_index(x, y)] << (x == nx - 1 ? "\n" : ",");
    }
  }
  ux_file.close();

  // Save uy_out
  std::ofstream uy_file(uy_filename);
  for (unsigned int y = 0; y < ny; ++y) {
    for (unsigned int x = 0; x < nx; ++x) {
      uy_file << uy_out[scalar_index(x, y)] << (x == nx - 1 ? "\n" : ",");
    }
  }
  uy_file.close();

  // Save rho_out
  std::ofstream rho_file(rho_filename);
  for (unsigned int y = 0; y < ny; ++y) {
    for (unsigned int x = 0; x < nx; ++x) {
      rho_file << rho_out[scalar_index(x, y)] << (x == nx - 1 ? "\n" : ",");
    }
  }
  rho_file.close();

  // Create XDMF file
  std::string xdmf_filename = "output_results/output_" + std::to_string(curr_time_step) + ".xdmf";
  std::ofstream xdmf_file(xdmf_filename);
  xdmf_file << "<?xml version=\"1.0\" ?>\n";
  xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xdmf_file << "<Xdmf Version=\"2.0\">\n";
  xdmf_file << "  <Domain>\n";
  xdmf_file << "    <Grid Name=\"LBM Grid\" GridType=\"Uniform\">\n";
  xdmf_file << "      <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << ny << " " << nx << "\"/>\n";
  xdmf_file << "      <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
  xdmf_file << "        <DataItem Dimensions=\"2\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n";
  xdmf_file << "          0.0 0.0\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "        <DataItem Dimensions=\"2\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n";
  xdmf_file << "          1.0 1.0\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Geometry>\n";

  // Write ux_out
  xdmf_file << "      <Attribute Name=\"ux_out\" AttributeType=\"Scalar\" Center=\"Node\">\n";
  xdmf_file << "        <DataItem Dimensions=\"" << ny << " " << nx << "\" Format=\"HDF\">\n";
  xdmf_file << "          " << ux_filename << "\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Attribute>\n";

  // Write uy_out
  xdmf_file << "      <Attribute Name=\"uy_out\" AttributeType=\"Scalar\" Center=\"Node\">\n";
  xdmf_file << "        <DataItem Dimensions=\"" << ny << " " << nx << "\" Format=\"HDF\">\n";
  xdmf_file << "          " << uy_filename << "\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Attribute>\n";

  // Write rho_out
  xdmf_file << "      <Attribute Name=\"rho_out\" AttributeType=\"Scalar\" Center=\"Node\">\n";
  xdmf_file << "        <DataItem Dimensions=\"" << ny << " " << nx << "\" Format=\"HDF\">\n";
  xdmf_file << "          " << rho_filename << "\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Attribute>\n";

  xdmf_file << "    </Grid>\n";
  xdmf_file << "  </Domain>\n";
  xdmf_file << "</Xdmf>\n";
  xdmf_file.close();
}
