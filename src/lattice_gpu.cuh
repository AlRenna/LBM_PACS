/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#ifndef __LATTICE_GPU_CUH__
#define __LATTICE_GPU_CUH__

#include "src/node.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

namespace lbm_gpu
{
  void cuda_simulation(unsigned int nx, 
                      unsigned int ny,
                      std::vector<Node> &nodes, // Changed from std::vector<NodeType> to std::vector<Node>
                      double tau,
                      double dt,
                      unsigned int save_iter,
                      unsigned int max_iter);
} //

#endif // __LATTICE_GPU_CUH__