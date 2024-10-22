/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

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


/**
 * @brief Class to model the lattice for the resolution of the LB method.
 * 
 */
class Lattice:
{

  public:

  /**
   * @brief Constructor
   * 
   * @param nx_ Lattice size (x)
   * @param ny_ Latice size (y)
   * @param nu_ Viscosity of the fluid
   */
  Lattice(unsigned int nx_, unsigned int ny_,
          double nu_);

  /**
   * @brief Function to set the initial and boundary conditions.
   * The function must be called from the main function before running the simulation so 
   * that the lattice can be properly initialized.
   * 
   * @param ux_in_ 
   * @param uy_in_ 
   * @param rho_in_ 
   * @param node_types_ 
   * @param boundary_node_dir_ 
   */
  void set_ICs_&_BCs(std::vector<double> ux_in_, 
               std::vector<double> uy_in_, 
               std::vector<double> rho_in_,
               std::vector<NodeType> node_types_,
               std::vector<std::vector<bool>> boundary_node_dir_);

  /**
   * @brief Function to populate the lattice nodes with the ICs and BCs. 
   * Also initialize the distribution functions to equilibrium.
   */
  void populate_Nodes();

  /**
   * @brief Run the simulation by iterating over time and over the lattice nodes
   */
  void run();



  /// @name Utility functions
  /// @{

  /**
   * @brief Function to compute the scalar index from the 2D coordinates.
   * 
   * @param x 
   * @param y 
   * @return size_t 
   */
  inline size_t scalar_index(unsigned int x, unsigned int y)
    {
      return nx*y+x;
    }

  ///@}

  /// @name Getters.
  /// @{

  /// Get the lattice size (nx) 
  inline unsigned int get_nx() const { return nx; }
  
  /// Get the lattice size (ny)
  inline unsigned int get_ny() const { return ny; }
  
  /// Get the relaxaion time 
  inline double get_tau() const { return tau; }
  
  /**
   * @brief Get the node object given the 2D coordinates.
   * 
   * @param x 
   * @param y 
   * @return Node& 
   */
  inline Node& get_node(unsigned int x, unsigned int y) const { return nodes[scalar_index(x, y)]; }

  /// @}
  
  private:
    /// @name Lattice variables
    /// @{

    /// Lattice size (x)
    const unsigned int nx;
    /// Lattice size (y)
    const unsigned int ny;

    /// Nodes in the lattice
    std::vector<Node> nodes;
    /// Node types (fluid, solid, etc.)
    std::vector<NodeType> node_types;
    /// Direction in which the BCs are to be applied for each boundary node 
    std::vector<std::vector<bool>> boundary_node_dir;

    /// @}

    /// @name Input-Output: physical quantities
    /// @{

    /// Input velocity
    std:vector<double> ux_in;
    /// Input velocity
    std:vector<double> uy_in;
    /// Input density
    std:vector<double> rho_in;

    /// Output velocity
    std:vector<double> ux_out;
    /// Output velocity
    std:vector<double> uy_out;
    /// Output density
    std:vector<double> rho_out;
    /// @}

    /// @name Physical parameters
    /// @{
    
    /// Viscosity
    const double nu; //= 1.0/6.0; 
    // Relaxation time
    const double tau; //= 3.0*nu+0.5;
    /// @}

    /// @name Time parameters
    /// @{
    
    /// Time step
    unsigned double dt = 0.1;
    
    /// Final time
    unsigned int T_final = 1;
    
    /// Maximum number of time steps
    unsigned int max_timesteps = 100;
    /// @}
    
}



#endif // __LATTICE_HPP__