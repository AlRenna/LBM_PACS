/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "lattice_gpu.cuh"

__global__ void hello(void)
{
  printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
}

void
lbm_gpu::cuda_simulation()
{
  std::cout << "Hello from CUDA!" << std::endl;
  hello <<<1, 1>>> ();
  cudaDeviceSynchronize(); // Ensure kernel execution is completed
}