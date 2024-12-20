/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#ifndef __LATTICE_GPU_CUH__
#define __LATTICE_GPU_CUH__

#include "src/node.hpp"
#include <omp.h>

#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>

namespace lbm_gpu
{
  /**
   * @brief Function to simulate the Lattice Boltzmann Method on the GPU.
   * 
   * @see lattice::run
   * @see lattice::run_gpu
   * 
   */
  void cuda_simulation(unsigned int nx, 
                      unsigned int ny,
                      std::vector<Node> &nodes, // Changed from std::vector<NodeType> to std::vector<Node>
                      double tau,
                      double dt,
                      double Cx,
                      double Crho,
                      unsigned int save_iter,
                      unsigned int max_iter);
} //

#endif // __LATTICE_GPU_CUH__