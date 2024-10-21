#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__

#include "node.hpp"

#include <vector>
#include <memory>

/*
Nella main function vogliamo richiamare uno script python che data un immagine in input
segmentata, ci restituisca:
- il vettore di NodeTypes
- le direzioni di boundary nodes
- 
*/

class Lattice:
{
  public:

  Lattice(unsigned int nx_, unsigned int ny_,
          double nu_);

  void set_ICs_&_BCs(std::vector<double> ux_in_, 
               std::vector<double> uy_in_, 
               std::vector<double> rho_in_,
               std::vector<NodeType> node_types_,
               std::vector<std::vector<bool>> boundary_node_dir_);

  void populate_Nodes();

  void run();



  // Utility functions

  inline size_t scalar_index(unsigned int x, unsigned int y)
    {
      return nx*y+x;
    }


  // Getter

  inline unsigned int get_nx() const { return nx; }
  inline unsigned int get_ny() const { return ny; }
  inline double get_tau() const { return tau; }
  inline Node& get_node(unsigned int x, unsigned int y) const { return nodes[y*nx+x]; }

  private:
    // Lattice size
    const unsigned int nx;
    const unsigned int ny;

    // Nodes
    std::vector<Node> nodes;
    std::vector<NodeType> node_types;
    std::vector<std::vector<bool>> boundary_node_dir;

    // Input-Output: physical quantities
    std:vector<double> ux_in;
    std:vector<double> uy_in;
    std:vector<double> rho_in;

    std:vector<double> ux_out;
    std:vector<double> uy_out;
    std:vector<double> rho_out;
    
    // Viscosity
    const double nu; //= 1.0/6.0; 
    // Relaxation time
    const double tau; //= 3.0*nu+0.5;
    
    // Time
    unsigned double dt = 0.1;
    unsigned int T_final = 1;
    unsigned int max_timesteps = 100;

    
}



#endif // __LATTICE_HPP__