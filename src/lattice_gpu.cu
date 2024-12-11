/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "src/lattice_gpu.cuh"

// __global__ void hello(void)
// {
//   printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
// }



void
lbm_gpu::cuda_simulation(unsigned int nx, 
                        unsigned int ny, 
                        std::vector<Node> &nodes,
                        double tau,
                        double dt,
                        unsigned int save_iter,
                        unsigned int max_iter)
{
  const int n = nx * ny;

  const int dir = Node::dir;

  //TODO: controlla per drag e lift se conviene fare in maniera diversa (somma in un unica variabile)

  // Host variables
  double * host_f_pre, * host_f_post, * host_f_adj, * host_ux, * host_uy, * host_rho, * host_drag, * host_lift, * host_bounce_back_delta;
  int * host_coord;
  bool * host_bounce_back_dir, * host_obstacle;
  NodeType * host_node_types;


  // Device variables
  double * d_f_pre, * d_f_post, * d_f_adj, * d_ux, * d_uy, * d_rho, * d_drag, * d_lift, * d_bounce_back_delta;
  int * d_coord;
  bool * d_bounce_back_dir, * d_obstacle;
  NodeType * d_node_types;

  // Allocate memory on the host
  cudaMalloc((void **) &host_f_pre, n * dir * sizeof(double));
  cudaMalloc((void **) &host_f_post, n * dir * sizeof(double));
  cudaMalloc((void **) &host_f_adj, n * dir * sizeof(double));
  cudaMalloc((void **) &host_ux, n * sizeof(double));
  cudaMalloc((void **) &host_uy, n * sizeof(double));
  cudaMalloc((void **) &host_rho, n * sizeof(double));
  cudaMalloc((void **) &host_drag, n * sizeof(double));
  cudaMalloc((void **) &host_lift, n * sizeof(double));
  cudaMalloc((void **) &host_coord, n * 2 * sizeof(int));
  cudaMalloc((void **) &host_bounce_back_delta, n * dir * sizeof(double));
  cudaMalloc((void **) &host_bounce_back_dir, n * dir * sizeof(bool));
  cudaMalloc((void **) &host_obstacle, sizeof(bool));
  cudaMalloc((void **) &host_node_types, n * sizeof(NodeType));

  // Allocate memory on the device
  cudaMalloc((void **) &d_f_pre, n * dir * sizeof(double));
  cudaMalloc((void **) &d_f_post, n * dir * sizeof(double));
  cudaMalloc((void **) &d_f_adj, n * dir * sizeof(double));
  cudaMalloc((void **) &d_ux, n * sizeof(double));
  cudaMalloc((void **) &d_uy, n * sizeof(double));
  cudaMalloc((void **) &d_rho, n * sizeof(double));
  cudaMalloc((void **) &d_drag, n * sizeof(double));
  cudaMalloc((void **) &d_lift, n * sizeof(double));
  cudaMalloc((void **) &d_coord, n * 2 * sizeof(int));
  cudaMalloc((void **) &d_bounce_back_delta, n * dir * sizeof(double));
  cudaMalloc((void **) &d_bounce_back_dir, n * dir * sizeof(bool));
  cudaMalloc((void **) &d_obstacle, sizeof(bool));
  cudaMalloc((void **) &d_node_types, n * sizeof(NodeType));


  // Set host data
  for(unsigned int index = 0; index < n; index++)
  {
    std::vector<double> temp_f_pre = nodes[index].get_f_pre();
    std::vector<double> temp_f_post = nodes[index].get_f_post();
    std::vector<double> temp_f_adj = nodes[index].get_f_adj();
    std::vector<double> temp_bounce_back_delta = nodes[index].get_bounce_back_delta();
    std::vector<bool> temp_bounce_back_dir = nodes[index].get_bounce_back_dir();

    host_ux[index] = nodes[index].get_ux();
    host_uy[index] = nodes[index].get_uy();
    host_rho[index] = nodes[index].get_rho();
    host_drag[index] = 0.0;
    host_lift[index] = 0.0;
    host_node_types[index] = nodes[index].get_node_type();
    
    if(host_node_types[index] == NodeType::obstacle && !(*host_obstacle)) {
      *host_obstacle = true;
    }

    for(unsigned int i = 0; i < dir; i++)
    {
      host_f_pre[index * dir + i] = temp_f_pre[i];
      host_f_post[index * dir + i] = temp_f_post[i];
      host_f_adj[index * dir + i] = temp_f_adj[i];
      host_bounce_back_delta[index * dir + i] = temp_bounce_back_delta[i];
      host_bounce_back_dir[index * dir + i] = temp_bounce_back_dir[i];
    }

    for(unsigned int i = 0; i < 2; i++)
    {
      host_coord[index * 2 + i] = nodes[index].get_coord()[i];
    }
    
  }

  // Copy data from host to device
  cudaMemcpy(d_f_pre, host_f_pre, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_f_post, host_f_post, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_f_adj, host_f_adj, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ux, host_ux, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_uy, host_uy, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rho, host_rho, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_drag, host_drag, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_lift, host_lift, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_coord, host_coord, n * 2 * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bounce_back_delta, host_bounce_back_delta, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bounce_back_dir, host_bounce_back_dir, n * dir * sizeof(bool), cudaMemcpyHostToDevice);
  cudaMemcpy(d_obstacle, host_obstacle, sizeof(bool), cudaMemcpyHostToDevice);
  cudaMemcpy(d_node_types, host_node_types, n * sizeof(NodeType), cudaMemcpyHostToDevice);

  // Free host memory
  // TODO: free host memory


  // Run simulation
  std::cout << "Running simulation\n" << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  unsigned int iter = 0;
  double total_time = 0.0;
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

  
  std::vector<double> vec_ux(nx * ny), vec_uy(nx * ny), vec_rho(nx * ny);
  vec_ux = arrayToVector(host_ux, nx * ny);
  vec_uy = arrayToVector(host_uy, nx * ny);
  vec_rho = arrayToVector(host_rho, nx * ny);
  // Save the initial conditions
  writeResults(u_file, ux_file, uy_file, rho_file, vec_ux, vec_uy, vec_rho, nx, ny);
}