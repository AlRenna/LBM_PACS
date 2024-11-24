/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__

#include "src/node.hpp"
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <vector>
#include <string>
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
 * The class contains the nodes of the lattice and the methods to run the simulation.
 * 
 */
class Lattice
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

    Lattice();

    /**
     * @brief Function to set the initial and boundary conditions.
     * The function must be called from the main function before running the simulation so 
     * that the lattice can be properly initialized.
     * 
     * @param ux_in_ Input velocity (x dirextion)
     * @param uy_in_ Input velocity (y direction)
     * @param rho_in_ Input density
     * @param filename_nodes File containing the information on the nodes of the lattice
     * 
     */
    void load_ICs_and_BCs(const std::vector<double>& ux_in_, 
                          const std::vector<double>& uy_in_, 
                          const std::vector<double>& rho_in_,
                          const std::string& filename_nodes);

    void readNodesFromCSV(const std::string& filename_nodes);


    /**
     * @brief Function to populate the lattice nodes with the ICs and BCs. 
     * Also initialize the distribution functions to equilibrium.
     */
    void populate_Nodes();

    /**
     * @brief Run the simulation by iterating over time and over the lattice nodes
     */
    void run();

    /**
     * @brief Function to write the results of the simulation to a file .
     * 
     * @param 
     */
    void writeResults(std::ofstream &file_u, std::ofstream &file_ux, std::ofstream &file_uy, std::ofstream &file_rho);

    /// @name Utility functions
    /// @{

    /**
     * @brief Function to compute the scalar index from the 2D coordinates.
     * 
     * @param x 
     * @param y 
     * @return size_t 
     */
    inline size_t scalar_index(unsigned int x, unsigned int y) const
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

    /// Get the time_step
    inline double get_dt() const { return dt; }
    
    /**
     * @brief Get the (const) node object given the 2D coordinates.
     * 
     * @param x 
     * @param y 
     * @return Node& 
     */
    inline const Node& get_node(unsigned int x, unsigned int y) const { return nodes[scalar_index(x, y)]; }
    
    /**
     * @brief Get the (non const) node object given the 2D coordinates.
     * 
     * @param x 
     * @param y 
     * @return Node& 
     */
    inline Node & get_node(unsigned int x, unsigned int y) { return nodes[scalar_index(x, y)]; }
    /// @}
  
  private:
    /// @name Lattice variables
    /// @{

    /// Lattice size (x)
    unsigned int nx;
    /// Lattice size (y)
    unsigned int ny;

    /// Nodes in the lattice
    std::vector<Node> nodes;
    /// Node types (fluid, solid, etc.)
    std::vector<NodeType> node_types;
    /// Directions in which the BCs are to be applied for each boundary node 
    std::vector<std::vector<bool>> bounce_back_dir;
    /// Distances between the node and the wall along the direction of bounce_back_dir for each boundary node
    std::vector<std::vector<double>> bounce_back_delta;

    /// @}

    /// @name Input-Output: physical quantities
    /// @{

    /// Input velocity
    std::vector<double> ux_in;
    /// Input velocity
    std::vector<double> uy_in;
    /// Input density
    std::vector<double> rho_in;

    /// Output velocity
    std::vector<double> ux_out;
    /// Output velocity
    std::vector<double> uy_out;
    /// Output density
    std::vector<double> rho_out;

    /// Output lift
    // double lift;
    /// Output drag
    // double drag;
    /// @}

    /// @name Physical parameters
    /// @{

    /// Characteristic Length 
    double Length;
    /// Viscosity
    double nu; //= 1.0/6.0; 
    // Relaxation time
    double tau; //= 3.0*nu+0.5;
    /// @}

    /// @name Time parameters
    /// @{
    
    /// Time step
    double dt;
    
    /// Final time
    double T_final = 1.;

    /// Save iteration (# of iterations after which the results are saved)
    unsigned int save_iter = 1;
    
    /// Maximum number of time steps
    unsigned int max_iter;
    /// @}
    
};



#endif // __LATTICE_HPP__