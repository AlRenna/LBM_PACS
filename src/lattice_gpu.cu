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

__device__ double compute_equilibrium(const double *d_weights, const double *d_coeff,
                                      double rho, double ux, double uy, int i)
{
  double weight = d_weights[i];
  double cx = d_coeff[2 * i];
  double cy = d_coeff[2 * i + 1];
  double u_dot_c = ux * cx + uy * cy;
  double u_sq = ux * ux + uy * uy;
  return weight * rho * (1.0 + 3.0 * u_dot_c + 4.5 * u_dot_c * u_dot_c - 1.5 * u_sq);
}

__device__ int find_forward_index(int current_index, int nx, int ny, int i, const double *d_coeff)
{
  int x = current_index % nx;
  int y = current_index / nx;
  double cx = d_coeff[2 * i];
  double cy = d_coeff[2 * i + 1];
  int x_new = x + cx;
  int y_new = y + cy;
  return y_new * nx + x_new;
}

__device__ int find_backward_index(int current_index, int nx, int ny, int i, const double *d_coeff, const int *d_bb_indexes)
{
  int x = current_index % nx;
  int y = current_index / nx;
  double cx = d_coeff[2 * d_bb_indexes[i]];
  double cy = d_coeff[2 * d_bb_indexes[i] + 1];
  int x_new = x + cx;
  int y_new = y + cy;
  return y_new * nx + x_new;
}

__device__ bool check_backward(int index, int nx, int ny, int i, const double *d_coeff, const int *d_bb_indexes, NodeType *d_node_types)
{
  int backward_index = find_backward_index(index, nx, ny, i, d_coeff, d_bb_indexes);
  return d_node_types[backward_index] == NodeType::fluid || d_node_types[backward_index] == NodeType::boundary;
}

__device__ void apply_IBB(const int dir,const double *d_weights, const double *d_coeff, const int *d_bb_indexes, 
                          double *d_f_post, double *d_f_adj,
                          double *d_ux, double *d_uy, double *d_rho,
                          NodeType *d_node_types, bool * d_bounce_back_dir, double * d_bounce_back_delta,
                          int nx, int ny, int i, int index)
{
  int forward_index = find_forward_index(index, nx, ny, i, d_coeff);
  double ux_wall = d_ux[forward_index];
  double uy_wall = d_uy[forward_index];

  if(check_backward(index, nx, ny, i, d_coeff, d_bb_indexes, d_node_types))
  {
    double cx = d_coeff[2 * i];
    double cy = d_coeff[2 * i + 1];

    double f_adj_post_coll = d_f_adj[index * dir + i];
    d_f_adj[index * dir + d_bb_indexes[i]] = (2 * d_bounce_back_delta[index * dir + i] * d_f_post[index * dir + i] + 
                    (1 - 2 * d_bounce_back_delta[index * dir + i]) * f_adj_post_coll) * 
                    (d_bounce_back_delta[index * dir + i] < 0.5) +
                    (1. / (2 * d_bounce_back_delta[index * dir + i]) * d_f_post[index * dir + i] + 
                    ((2 * d_bounce_back_delta[index * dir + i] - 1.) / (2 * d_bounce_back_delta[index * dir + i])) * d_f_post[index * dir + d_bb_indexes[i]]) *
                    (d_bounce_back_delta[index * dir + i] >= 0.5) - 
                    (ux_wall * cx + uy_wall * cy) * d_weights[i] * 6;
  }
  else
  {
    d_f_adj[index * dir + d_bb_indexes[i]] = d_f_post[index * dir + i];
  }
}

__device__ void apply_anti_BB(const int dir,const double *d_weights, const double *d_coeff, const int *d_bb_indexes, 
                          double *d_f_post, double *d_f_adj,
                          double *d_ux, double *d_uy, double *d_rho,
                          NodeType *d_node_types, bool * d_bounce_back_dir, double * d_bounce_back_delta,
                          int nx, int ny, int i, int index)
{
  int forward_index = find_forward_index(index, nx, ny, i, d_coeff);
  double ux_wall = d_ux[forward_index];
  double uy_wall = d_uy[forward_index];
  double rho_wall = d_rho[forward_index];

  if(d_node_types[forward_index] == NodeType::outlet)
  {
    int backward_index = find_backward_index(index, nx, ny, i, d_coeff, d_bb_indexes);
    double ux_fluid = d_ux[backward_index];
    double uy_fluid = d_uy[backward_index];

    ux_wall = (d_ux[forward_index] + ux_fluid) / 2;
    uy_wall = (d_uy[forward_index] + uy_fluid) / 2;
    
    rho_wall = 0.8 * (2 * (d_f_post[index * dir + 1] + d_f_post[index * dir + 5] + d_f_post[index * dir + 8]) + d_f_post[index * dir + 0] + d_f_post[index * dir + 2] + d_f_post[index * dir + 4]) / (1. - ux_wall);
  }

  if(check_backward(index, nx, ny, i, d_coeff, d_bb_indexes, d_node_types))
  {
    double cx = d_coeff[2 * i];
    double cy = d_coeff[2 * i + 1];

    d_f_adj[index * dir + d_bb_indexes[i]] = -d_f_post[index * dir + i] +
                          2 * d_weights[i] * rho_wall *
                          (1 + 4.5 * (cx * ux_wall + cy * uy_wall) * (cx * ux_wall + cy * uy_wall) -
                          3.5 * (ux_wall * ux_wall + uy_wall * uy_wall));
  }
  else
  {
    d_f_adj[index * dir + d_bb_indexes[i]] = d_f_post[index * dir + i];
  }
}

__global__ void collide_and_stream_kernel(
  const int dir, const double *d_weights, const double *d_coeff, const int *d_bb_indexes,
  double *d_f_pre, double *d_f_post, double *d_f_adj, 
  double *d_ux, double *d_uy, double *d_rho, 
  NodeType *d_node_types, bool * d_bounce_back_dir, int nx, int ny, double tau) 
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int n = nx * ny;

  if(index < n) 
  {
    if(d_node_types[index] == NodeType::fluid || d_node_types[index] == NodeType::boundary) 
    {
      // Collision step
      for(int i = 0; i < dir; ++i)
      {
        double feq = compute_equilibrium(d_weights, d_coeff, d_rho[index], d_ux[index], d_uy[index], i);
        d_f_post[index * dir + i] = d_f_pre[index * dir + i] - (d_f_pre[index * dir + i] - feq) / tau;
      }

      // Streaming step
      d_f_adj[index * dir + 0] = d_f_post[index * dir + 0];
      for(int i = 0; i < dir; ++i) 
      {
        if(!d_bounce_back_dir[index * dir + i]) {
          int index_new = find_forward_index(index, nx, ny, i, d_coeff);
          d_f_adj[index_new * dir + i] = d_f_post[index * dir + i];
        }
      }
    }
  }
}

__global__ void apply_BCs_and_compute_quantities_kernel(
  const int dir, const double *d_weights, const double *d_coeff, const int *d_bb_indexes,
  double *d_f_pre, double *d_f_post, double *d_f_adj,
  double *d_ux, double *d_uy, double *d_rho,
  double *d_drag, double *d_lift, bool * d_obstacle,
  NodeType *d_node_types, bool * d_bounce_back_dir, double * d_bounce_back_delta, int nx, int ny) 
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int n = nx * ny;

  if(index < n)
  {
    if(d_node_types[index] == NodeType::fluid || d_node_types[index] == NodeType::boundary)
    {
      if(d_node_types[index] == NodeType::boundary)
      {
        // Apply boundary conditions
        for(int i = 0; i < dir; ++i)
        {
          if(d_bounce_back_dir[index * dir + i])
          {
            int index_new = find_forward_index(index, nx, ny, i, d_coeff);

            if(d_node_types[index_new] == NodeType::solid ||
              d_node_types[index_new] == NodeType::obstacle ||
              d_node_types[index_new] == NodeType::inlet)
            {
              // Interpolated Bounce-Back
              apply_IBB(dir, d_weights, d_coeff, d_bb_indexes, d_f_post, d_f_adj, d_ux, d_uy, d_rho, d_node_types, d_bounce_back_dir, d_bounce_back_delta, nx, ny, i, index);
            }
            else if(d_node_types[index_new] == NodeType::outlet)
            {
              // Anti Bounce-Back
              apply_anti_BB(dir, d_weights, d_coeff, d_bb_indexes, d_f_post, d_f_adj, d_ux, d_uy, d_rho, d_node_types, d_bounce_back_dir, d_bounce_back_delta, nx, ny, i, index);
            }
            else
            {
              printf("Error: Invalid BCs type at index %d\n", index);
              return;
            }
          }
        }

        // Compute drag and lift
        if(*d_obstacle)
        {
          double dr = 0.0;
          double lf = 0.0;
          for(int i = 0; i < dir; ++i)
          {
            int forward_index = find_forward_index(index, nx, ny, i, d_coeff);
            if(d_node_types[forward_index] == NodeType::obstacle)
            {
              double cx = d_coeff[2 * i];
              double cy = d_coeff[2 * i + 1];
              double cx_bb = d_coeff[2 * d_bb_indexes[i]];
              double cy_bb = d_coeff[2 * d_bb_indexes[i] + 1];

              dr += cx * d_f_pre[index * dir + i] - cx_bb * d_f_adj[index * dir + d_bb_indexes[i]];
              lf += cy * d_f_pre[index * dir + i] - cy_bb * d_f_adj[index * dir + d_bb_indexes[i]];
            }
          }
          atomicAdd(d_drag, dr);
          atomicAdd(d_lift, lf);
        }
      }

      // Update f
      for(int i = 0; i < dir; ++i)
      {
        d_f_pre[index * dir + i] = d_f_adj[index * dir + i];
      }

      // Compute macroscopic quantities
      double rho = 0.0;
      double ux = 0.0;
      double uy = 0.0;
      for(int i = 0; i < dir; ++i)
      {
        double f = d_f_pre[index * dir + i];
        rho += f;
        ux += f * d_coeff[2 * i];
        uy += f * d_coeff[2 * i + 1];
      }
      d_rho[index] = rho;
      d_ux[index] = ux / rho;
      d_uy[index] = uy / rho;
    }
  }
}

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

  // Constants for CUDA kernel
  const int dir = Node::dir;
  const double *weights = vectorToArray(Node::weights);
  const double *coeff = vector2DToArray(Node::coeff);
  const int *bb_indexes = vectorToArray(Node::bb_indexes);

  // Copy constant variables to device
  double *d_weights, *d_coeff;
  int *d_bb_indexes;

  cudaMalloc((void **) &d_weights, dir * sizeof(double));
  cudaMalloc((void **) &d_coeff, 2 * dir * sizeof(double));
  cudaMalloc((void **) &d_bb_indexes, dir * sizeof(int));

  cudaMemcpy(d_weights, weights, dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_coeff, coeff, 2 * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bb_indexes, bb_indexes, dir * sizeof(int), cudaMemcpyHostToDevice);


  std::vector<double> lift_out(max_iter, 0.0);
  std::vector<double> drag_out(max_iter, 0.0);
  
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
  cudaMallocHost((void **) &host_f_pre, n * dir * sizeof(double));
  cudaMallocHost((void **) &host_f_post, n * dir * sizeof(double));
  cudaMallocHost((void **) &host_f_adj, n * dir * sizeof(double));
  cudaMallocHost((void **) &host_ux, n * sizeof(double));
  cudaMallocHost((void **) &host_uy, n * sizeof(double));
  cudaMallocHost((void **) &host_rho, n * sizeof(double));
  cudaMallocHost((void **) &host_drag, sizeof(double));
  cudaMallocHost((void **) &host_lift, sizeof(double));
  cudaMallocHost((void **) &host_coord, n * 2 * sizeof(int));
  cudaMallocHost((void **) &host_bounce_back_delta, n * dir * sizeof(double));
  cudaMallocHost((void **) &host_bounce_back_dir, n * dir * sizeof(bool));
  cudaMallocHost((void **) &host_obstacle, sizeof(bool));
  cudaMallocHost((void **) &host_node_types, n * sizeof(NodeType));

  // Initialize host_obstacle
  *host_obstacle = false;

  // Allocate memory on the device
  cudaMalloc((void **) &d_f_pre, n * dir * sizeof(double));
  cudaMalloc((void **) &d_f_post, n * dir * sizeof(double));
  cudaMalloc((void **) &d_f_adj, n * dir * sizeof(double));
  cudaMalloc((void **) &d_ux, n * sizeof(double));
  cudaMalloc((void **) &d_uy, n * sizeof(double));
  cudaMalloc((void **) &d_rho, n * sizeof(double));
  cudaMalloc((void **) &d_drag, sizeof(double));
  cudaMalloc((void **) &d_lift, sizeof(double));
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
  *host_drag = 0.0;
  *host_lift = 0.0;

  std::cout << "Copying data to device\n" << std::endl;
  // Copy data from host to device
  cudaMemcpy(d_f_pre, host_f_pre, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_f_post, host_f_post, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_f_adj, host_f_adj, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ux, host_ux, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_uy, host_uy, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rho, host_rho, n * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_drag, host_drag, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_lift, host_lift, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_coord, host_coord, n * 2 * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bounce_back_delta, host_bounce_back_delta, n * dir * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_bounce_back_dir, host_bounce_back_dir, n * dir * sizeof(bool), cudaMemcpyHostToDevice);
  cudaMemcpy(d_obstacle, host_obstacle, sizeof(bool), cudaMemcpyHostToDevice);
  cudaMemcpy(d_node_types, host_node_types, n * sizeof(NodeType), cudaMemcpyHostToDevice);

  // Free host memory
  cudaFreeHost(host_f_pre);
  cudaFreeHost(host_f_post);
  cudaFreeHost(host_f_adj);
  cudaFreeHost(host_coord);
  cudaFreeHost(host_bounce_back_delta);
  cudaFreeHost(host_bounce_back_dir);
  cudaFreeHost(host_obstacle);
  cudaFreeHost(host_node_types);


  // Run simulation
  std::cout << "Running simulation\n" << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  unsigned int iter = 0;
  double total_time = 0.0;
  std::cout << "Create folder and files\n" << std::endl;
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

  std::cout << "Save initial conditions\n" << std::endl;
  std::vector<double> vec_ux(nx * ny), vec_uy(nx * ny), vec_rho(nx * ny);
  vec_ux = arrayToVector(host_ux, nx * ny);
  vec_uy = arrayToVector(host_uy, nx * ny);
  vec_rho = arrayToVector(host_rho, nx * ny);
  // Save the initial conditions
  writeResults(u_file, ux_file, uy_file, rho_file, vec_ux, vec_uy, vec_rho, nx, ny);

  iter = iter + 1;

  std::cout << "Start simulation loop\n" << std::endl;
  while(iter <= max_iter) {
    auto iter_start_time = std::chrono::high_resolution_clock::now();

    // if(iter % save_iter == 0 || iter == max_iter - 1) {
      std::cout << "Iteration: " << iter << std::endl;
      std::cout << "Time: " << iter * dt << std::endl;
      std::cout << "Collision and streaming" << std::endl;
    // }

    // Define block size
    int blockSize = 256; // 256 or 512

    // Calculate grid size
    int gridSize = (nx * ny + blockSize - 1) / blockSize;

    // Launch CUDA kernel for collision and streaming
    collide_and_stream_kernel<<<gridSize, blockSize>>>(dir, d_weights, d_coeff, d_bb_indexes, 
                                                      d_f_pre, d_f_post, d_f_adj,
                                                      d_ux, d_uy, d_rho, 
                                                      d_node_types, d_bounce_back_dir, nx, ny, tau);

    // if(iter % save_iter == 0 || iter == max_iter - 1) {
      std::cout << "Physical quantities evaluation\n" << std::endl;
    // }

    // Launch CUDA kernel for applying boundary conditions and computing physical quantities
    apply_BCs_and_compute_quantities_kernel<<<gridSize, blockSize>>>(dir, d_weights, d_coeff, d_bb_indexes, 
                                                                    d_f_pre, d_f_post, d_f_adj,
                                                                    d_ux, d_uy, d_rho,
                                                                    d_drag, d_lift, d_obstacle,
                                                                    d_node_types, d_bounce_back_dir, d_bounce_back_delta, nx, ny);

    // Copy lift and drag results from device to host
    cudaMemcpy(&host_lift, d_lift, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_drag, d_drag, sizeof(double), cudaMemcpyDeviceToHost);

    lift_out[iter] = *host_lift;
    drag_out[iter] = *host_drag;

    if(iter % save_iter == 0 || iter == max_iter - 1) {
      // Copy results from device to host
      cudaMemcpy(host_ux, d_ux, n * sizeof(double), cudaMemcpyDeviceToHost);
      cudaMemcpy(host_uy, d_uy, n * sizeof(double), cudaMemcpyDeviceToHost);
      cudaMemcpy(host_rho, d_rho, n * sizeof(double), cudaMemcpyDeviceToHost);

      vec_ux = arrayToVector(host_ux, nx * ny);
      vec_uy = arrayToVector(host_uy, nx * ny);
      vec_rho = arrayToVector(host_rho, nx * ny);

      writeResults(u_file, ux_file, uy_file, rho_file, vec_ux, vec_uy, vec_rho, nx, ny);
    }

    iter = iter + 1;
    auto iter_end_time = std::chrono::high_resolution_clock::now();
    total_time += std::chrono::duration<double>(iter_end_time - iter_start_time).count();
  }

  // Save the lift and drag results
  if(host_obstacle)
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
  std::chrono::duration<double> elapsed_time = end_time - start_time;
  double mean_time_per_iter = total_time / max_iter;
  std::cout << "Simulation completed in " << elapsed_time.count() << " seconds" << std::endl;
  std::cout << "Mean time per iteration: " << mean_time_per_iter << " seconds.\n" << std::endl;

  u_file.close();
  ux_file.close();
  uy_file.close();
  rho_file.close();

  // Free device memory
  cudaFreeHost(host_ux);
  cudaFreeHost(host_uy);
  cudaFreeHost(host_rho);
  cudaFreeHost(host_drag);
  cudaFreeHost(host_lift);
}