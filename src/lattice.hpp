/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__

#include "src/node.hpp"
#include "src/utils/utils.hpp"

#include <nlohmann/json.hpp>
#include <omp.h>
#include <chrono>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <vector>
#include <string>
#include <memory>


/**
 * @brief Class to model the lattice for the resolution of the LB method. 
 * The class contains the nodes of the lattice and the methods to run the simulation.
 * 
 */
class Lattice
{
  public:

    /**
     * @brief Constructor. Reads the parameters from params.json file.
     * 
     */
    Lattice();

    /**
     * @brief Function to initialize the lattice.
     * It sets the input velocity, density, reads the lattice information from a .csv file, 
     * and populates the lattice nodes.
     * The function must be called from the main function before running the simulation so 
     * that the lattice can be properly initialized.
     * 
     * @param ux_in_ Input velocity (x dirextion)
     * @param uy_in_ Input velocity (y direction)
     * @param rho_in_ Input density
     * @param filename_nodes File containing the information on the nodes of the lattice
     * 
     */
    void initialize(const std::vector<double>& ux_in_, 
                          const std::vector<double>& uy_in_, 
                          const std::vector<double>& rho_in_,
                          const std::string& filename_nodes);
    
    /**
     * @brief Function to popolate the member of the lattice by reading a .cvs file containing the lattice information.
     * 
     * @param filename_nodes 
     */
    void readNodesFromCSV(const std::string& filename_nodes);

    /**
     * @brief Function to populate the lattice nodes with the ICs and BCs. 
     * Also initialize the distribution functions to equilibrium.
     */
    void populate_Nodes();

    /**
     * @brief Call the run (with parallelization on CPU or GPU) function to start the simulation.
     * 
     * If the first argument is "-gpu" the simulation is run on the GPU, otherwise on the CPU.
     */
    void run(int argc, char **argv);

    /**
     * @brief Run the simulation (with OpenMP parallelization)
     *  
     */
    void run_cpu();

    /**
     * @brief Run the simulation on the GPU (with CUDA parallelization)
     * 
     */
    void run_gpu();

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

    /**
     * @brief Function to compute the 2D coordinates from the scalar index.
     * 
     * @param index 
     * @return std::pair<unsigned int, unsigned int> 
     */
    inline std::vector<unsigned int> coord_index(size_t index) const
      {
        return {static_cast<unsigned int>(index % nx), static_cast<unsigned int>(index / nx)};
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

    /// Get the inverse relaxation time
    inline double get_tauinv() const { return tauinv; }

    /// Get the 1 - inverse relaxation time
    inline double get_omtauinv() const { return omtauinv; }

    /// Get the time_step
    inline double get_dt() const { return dt; }

    /// Get simulate time
    inline double get_time() const { return static_cast<double>(iter) / max_iter; }
    
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
    std::vector<double> lift_out;
    /// Output drag
    std::vector<double> drag_out;
    /// @}

    /// @name Physical parameters
    /// @{

    /// Viscosity
    double nu;
    /// Relaxation time
    double tau;
    /// Inverse relaxation time
    double tauinv;
    /// 1 - inverse relaxation time
    double omtauinv;

    /// @}

    /// @name Time parameters
    /// @{
    
    /// Time step
    double dt;
    
    /// Final time
    double T_final = 1.;

    /// Save time
    double save_time = 0.1;

    /// Save iteration (# of iterations after which the results are saved)
    unsigned int save_iter = 1;

    /// iter
    unsigned int iter;
    
    /// Maximum number of time steps
    unsigned int max_iter;
    /// @}

    /// @name Conversion variables 
    /// @{

    /// Conversion factor from lattice space units to physical units
    double Cx;

    /// Conversion factor from lattice density units to physical units
    double Crho;

    /// @}
    
};



#endif // __LATTICE_HPP__