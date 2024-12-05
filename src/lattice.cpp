/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "src/lattice.hpp"

//TODO: Fare un esempio con ala (ostacolo)

Lattice::Lattice()
{
  std::ifstream param_file("params.json");
  if (!param_file.is_open()) {
    throw std::runtime_error("Could not open params.json");
  }

  nlohmann::json param_json;
  param_file >> param_json;
  param_file.close();

  if (param_json.find("lattice") == param_json.end() || 
      param_json.find("generated_variables") == param_json.end() ||
      param_json["generated_variables"].find("new_nx") == param_json["generated_variables"].end() || 
      param_json["generated_variables"].find("new_ny") == param_json["generated_variables"].end() || 
      param_json["generated_variables"].find("iterations") == param_json["generated_variables"].end() ||
      param_json["lattice"].find("nu") == param_json["lattice"].end()) {
    throw std::runtime_error("param.json does not contain required parameters");
  }

  Length = param_json["lattice"]["Length"];
  nx = param_json["generated_variables"]["new_nx"];
  ny = param_json["generated_variables"]["new_ny"];
  dt = param_json["generated_variables"]["dt"];
  max_iter = param_json["generated_variables"]["iterations"];
  nu = param_json["lattice"]["nu"];  
  T_final = param_json["time"]["T_final"];
  save_time = param_json["time"]["save_time"];

  save_iter = static_cast<unsigned int>(std::ceil(save_time / dt));
  
  // Cs = 1/sqrt(3) * delta_x / delta_t = dx/dt; 
  // dx = L/n; n = sqrt(nx^2 + ny^2); dt = sqrt(Cs) * dx;
  // max_iter = T_final / dt; 
  // dt = std::sqrt(3) * Length / std::sqrt(nx * nx + ny * ny);
  // max_iter = static_cast<int>(std::ceil(T_final / dt));

  tau = 3.0 * nu + 0.5;
  nodes.resize(nx*ny);
  node_types.resize(nx*ny, NodeType::solid);
  bounce_back_delta.resize(nx*ny, std::vector<double>(Node::dir, 0.0));
  bounce_back_dir.resize(nx*ny, std::vector<bool>(Node::dir, false));
  ux_in.resize(nx*ny, 0.);
  uy_in.resize(nx*ny, 0.);
  rho_in.resize(nx*ny, 1.0);
  ux_out.resize(nx*ny, 0.);
  uy_out.resize(nx*ny, 0.);
  rho_out.resize(nx*ny, 1.0);
  lift_out.resize(max_iter, 0.);
  drag_out.resize(max_iter, 0.);
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

  // Find Boundary Nodes types

  // Populate the nodes
  populate_Nodes();
}

void 
Lattice::readNodesFromCSV(const std::string& filename) 
{
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
    
    if(nodeType == NodeType::boundary){
      // Read the last 8 columns
      for (int i = 1; i < Node::dir; ++i) {
        std::getline(ss, token, ',');
        double value = std::stod(token);
        bounce_back_delta[index][i] = value;
        bounce_back_dir[index][i] = (value > std::numeric_limits<double>::epsilon());
      }
    }
  }

  file.close();
}

void
Lattice::populate_Nodes()
{
  std::cout << "Populating Nodes" << std::endl;

  #pragma omp parallel for
  for(unsigned int index = 0; index<nx*ny; ++index){
    std::vector<unsigned int> coord_xy = coord_index(index);
    nodes[index] = Node(node_types[index], coord_xy, ux_in[index], uy_in[index], rho_in[index]);
    if(node_types[index] == NodeType::boundary){
      nodes[index].set_bounce_back_properties(bounce_back_dir[index], bounce_back_delta[index]);
    }
    nodes[index].init_equilibrium();
  }
}

void
Lattice::run(int argc, char **argv)
{
  bool gpuFlag = false;
  if (argc > 1 && std::string(argv[1]) == "-gpu") {
    gpuFlag = true;
  }

  #ifdef USE_CUDA
    if (gpuFlag)
    { 
        std::cout << "CUDA enabled. Running on GPU" << std::endl;
        run_gpu();
    }
    else
    {   
        std::cout << "CUDA enabled. Running on CPU" << std::endl;
        run_cpu();
    }
  #else
    std::cout << "CUDA not enabled. Running on CPU" << std::endl;
    //run_cpu();
  #endif
}

void
Lattice::run_cpu()
{
  std::cout << "Running simulation\n" << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  unsigned int iter = 0;
  double total_time = 0.0;

  // Check if there is an obstacle
  bool obstacle_present = false;
  for(unsigned int index = 0; index<nx*ny; ++index){
    if(node_types[index] == NodeType::obstacle){
      obstacle_present = true;
      std::cout << "Obstacle present\n" << std::endl;
      break;
    }
  }

  // Delete the output_results directory if it exists
  if (std::filesystem::exists("output_results")) {
    std::filesystem::remove_all("output_results");
    std::filesystem::create_directory("output_results");
  }
  else{
    std::filesystem::create_directory("output_results");
  }

  if (std::filesystem::exists("output_animations")) {
    std::filesystem::remove_all("output_animations");
    std::filesystem::create_directory("output_animations");
  }
  else{
    std::filesystem::create_directory("output_animations");
  }
  
  std::string u_filename = "output_results/velocity_out.txt";
  std::string ux_filename = "output_results/ux_out.txt";
  std::string uy_filename = "output_results/uy_out.txt";
  std::string rho_filename = "output_results/rho_out.txt";

  std::ofstream u_file(u_filename);
  std::ofstream ux_file(ux_filename);
  std::ofstream uy_file(uy_filename);
  std::ofstream rho_file(rho_filename);
  

  // Save the initial conditions
  writeResults(u_file, ux_file, uy_file, rho_file);
  //writeResults(iter);

  iter = iter + 1;
  
  while(iter <= max_iter)
  {
    auto iter_start_time = std::chrono::high_resolution_clock::now();

    if( iter%save_iter == 0 || iter == max_iter-1)
    {
      std::cout << "Iteration: " << iter << std::endl;
      std::cout << "Time: " << iter * dt << std::endl;
      std::cout << "Collision and streaming" << std::endl;
    }

    #pragma omp parallel for 
    for(unsigned int index = 0; index<nx*ny; ++index){
      if(node_types[index] == NodeType::fluid || node_types[index] == NodeType::boundary)
      {
        nodes[index].collide(*this);
        nodes[index].stream(*this);
      }
    }

    if( iter%save_iter == 0 || iter == max_iter-1) 
      std::cout << "Physical quantities evaluation\n" << std::endl;

    double lift = 0.0;
    double drag = 0.0;
    #pragma omp parallel for default(none) \
            reduction(+:lift) reduction(+:drag) shared(obstacle_present)
    for(unsigned int index = 0; index<nx*ny; ++index){
      if(node_types[index] == NodeType::fluid || node_types[index] == NodeType::boundary)
      {
        if(node_types[index] == NodeType::boundary){
          nodes[index].apply_BCs(*this);
          if(obstacle_present)
          { 
            nodes[index].compute_drag_and_lift(*this);
            lift = lift + nodes[index].get_lift();
            drag = drag + nodes[index].get_drag();  
          }
        }
        nodes[index].update_f();

        nodes[index].compute_physical_quantities();

        ux_out[index] = nodes[index].get_ux();
        uy_out[index] = nodes[index].get_uy();
        rho_out[index] = nodes[index].get_rho();
                
      }
    }
    lift_out[iter] = lift;
    drag_out[iter] = drag;

    if( iter%save_iter == 0 || iter == max_iter-1)
    {
      writeResults(u_file, ux_file, uy_file, rho_file);
    }
    
    iter = iter + 1;
    auto iter_end_time = std::chrono::high_resolution_clock::now();
    total_time += std::chrono::duration<double>(iter_end_time - iter_start_time).count();
  }

  // Save the lift and drag results
  if(obstacle_present)
  { 
    std::string lift_drag_filename = "output_results/lift_&_drag.txt";
    std::ofstream lift_drag_file(lift_drag_filename);
    lift_drag_file << "Lift:\n";
    // Save the lift and drag
    for(unsigned int t=0; t<max_iter; ++t){
      lift_drag_file << lift_out[t] << " ";
    }
    
    lift_drag_file << "\nDrag:\n";
    // Save the lift and drag
    for(unsigned int t=0; t<max_iter; ++t){
      lift_drag_file << drag_out[t] << " ";
    }
    lift_drag_file.close();
  }


  auto end_time = std::chrono::high_resolution_clock::now();
  double elapsed_time = std::chrono::duration<double>(end_time - start_time).count();
  double mean_time_per_iter = total_time / max_iter;
  std::cout << "Simulation completed in " << elapsed_time << " seconds." << std::endl;
  std::cout << "Mean time per iteration: " << mean_time_per_iter << " seconds.\n" << std::endl;

  u_file.close();
  ux_file.close();
  uy_file.close();
  rho_file.close();
  
}


#ifdef USE_CUDA
#include "src/lattice_gpu.cuh"
void
Lattice::run_gpu()
{
  lbm_gpu::cuda_simulation();
}
#endif

void 
Lattice::writeResults(std::ofstream &file_u, std::ofstream &file_ux, std::ofstream &file_uy, std::ofstream &file_rho) {
  // Save ux_out
  for(unsigned int index = 0; index<nx*ny; ++index){
    file_u << std::sqrt(ux_out[index] * ux_out[index] + uy_out[index] * uy_out[index]) << " ";
    file_ux<< ux_out[index] << " ";
    file_uy<< uy_out[index] << " ";
    file_rho << rho_out[index] << " ";
  }
  file_u << "\n";
  file_ux << "\n";
  file_uy << "\n";
  file_rho << "\n";
  
}
