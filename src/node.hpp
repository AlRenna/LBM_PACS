/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#ifndef __NODE_HPP__
#define __NODE_HPP__

#include <vector>
#include <memory>
#include <cmath>

#include "utils.cpp"

// Forward declaration
class Lattice;

/**
 * @brief Enumeration class for the node types. Used to handle the different Initialization and Boundary conditions.
 * 
 */
enum class NodeType
{
  fluid,
  solid,
  inflow,
  outflow,
  boundary
};

/**
 * @brief Node class. Contains the information of the node in the lattice.
 * 
 */
class Node
{
  public:
    /**
     * @brief Constructor
     * 
     * @param node_type_ 
     * @param coord_ 
     * @param ux_ 
     * @param uy_ 
     * @param rho_ 
     */
    Node(NodeType node_type_, std::vector<unsigned int> coord_,
        double ux_, double uy_, double rho_);

    Node() = default;
    // ~Node();

    
    
    /**
     * @brief Initialize the distribution functions to the equilibrium values.
     */
    void init_equilibrium();

    /**
     * @brief Collision step.
     * Compute and relax to the equilibrium distribution. The BGK collision operator is used: 
     * f_i(t+dt) = f_i(t) - 1/tau (f_i(t) - f_i^eq(t))
     * 
     * @param lattice 
     */
    void collide(const Lattice &lattice);

    /**
     * @brief Streaming step.
     * Stream the distribution functions to the adjacent nodes.
     * 
     */
    void stream(Lattice& lattice);

     /**
     * @brief Apply the Iterpolated Bounce-Back.
     * Use the post-collision distribution to compute the IBB:
     * 
     * f_bb_i = 2*q_i*f_i(x_b) + (1 - 2*q)*f_i(x_f) if q_i < 1/2
     * 
     * f_bb_i = 1/(2*q_i)*f_i(x_b) + [(2*q - 1)/(2*q)]*f_bb_i(x_b) if q_i >= 1/2
     *
     */
    void apply_IBB();


    /**
     * @brief Compute the physical quantities (velocity components and density) from the distribution functions.
     */
    void compute_physical_quantities();


    void update_f();



    /// @name Getters
    /// @{ 
    
    /**
     * @brief Get the distribution function at index (direction) i.
     */
    inline double get_f(int i) const { return (*f)[i]; }
    
    /// @}

    /// @name Setters
    /// @{

    /**
     * @brief Fill the vectors containg the boundary node properties, 
     * such as the direction in which to apply the BCs (IBB) and the 
     * weigthed distance between the node and the wall.
     * 
     * @param boundary_node_dir_ 
     */
    void set_boundary_node_properties(std::vector<bool> boundary_node_dir_, 
                                      std::vector<double> boundary_node_delta_,
                                      double dt);

    /**
     * @brief Set the distribution for the adjacent node alorng direction i.
     * 
     * @param i 
     * @param value 
     */
    inline void set_f_adj(int i, double value) { (*f_adj)[i] = value; }
    /// @}



    /// @name Static lattice variables "D2Q9 model"
    /// @{

    /// Number of dimensions
    static const int dim  = 2;
    
    /**
     * @brief Number of velocity directions
     * Direction numbering scheme:
     * 
     * 6 2 5
     * 
     * 3 0 1
     * 
     * 7 4 8
     */
    static const int dir = 9;

    /// Weights
    static const std::vector<double> weights;
    /// Coefficients
    static const std::vector<std::vector<double>> coeff;
    
    /// Opposite direction indexes for Bounce Back
    static const std::vector<int> bb_indexes;

    /// @}

  private:

    /**
     * @brief Velocity distribution functions.
     * Direction numbering scheme:
     * 
     * 6 2 5
     * 
     * 3 0 1
     * 
     * 7 4 8
     */
    std::unique_ptr<std::vector<double>> f;

    /// Velocity distribution functions from adjacent nodes
    std::unique_ptr<std::vector<double>> f_adj;

    /// Node velocity components (x)
    double ux;
    /// Node velocity components (y)
    double uy;
    /// Node density
    double rho;

    /// Node type (fluid, boundary, solid)
    NodeType node_type;
    /// Directions in which the BCs are to be applied around the node (if it's a boundary node)
    std::vector<bool> boundary_node_dir;
    /// Distance between the node and the wall (if it's a boundary node), along the direction of boundary_node_dir
    std::vector<double> boundary_node_delta;
    /// 2D coordinates of the node
    std::vector<unsigned int> coord;


};




#endif // __NODE_HPP__