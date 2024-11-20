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

#include "src/utils/utils.hpp"

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
  inlet,
  outlet,
  boundary
};

/**
 * @brief Node class. Contains the information of the node in the lattice such as: type, position, velocity
 * distribution and physical quantities.
 * 
 */
class Node
{
  public:
    /**
     * @brief Constructor
     * 
     * @param node_type_ A NodeType enumerator indicating the type of the node (fluid, boundary, solid).
     * @param coord_ 2D coordinates of the node.
     * @param ux_ Initial velocity component in the x direction.
     * @param uy_ Initial velocity component in the y direction.
     * @param rho_ Initial density of the node.
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
    void apply_IBB(const Lattice &lattice);

    void apply_BB(const Lattice &lattice);


    /**
     * @brief Compute the physical quantities (velocity components and density) from the distribution functions.
     */
    void compute_physical_quantities();

    /**
     * @brief Function to compute physical integrals such as lift and drag
     * 
     */
    // void compute_integrals();


    void update_f();



    /// @name Getters
    /// @{ 
    
    /**
     * @brief Get the distribution function at index (direction) i.
     */
    inline double get_f(int i) const { return (*f)[i]; }

    inline double get_ux() const { return ux; }
    inline double get_uy() const { return uy; }
    inline double get_rho() const { return rho; }
    
    /// @}

    /// @name Setters
    /// @{

    /**
     * @brief Fill the vectors containg the boundary node properties, 
     * such as the direction in which to apply the BCs (IBB) and the 
     * weigthed distance between the node and the wall.
     * 
     * @param bounce_back_dir_ 
     */
    void set_bounce_back_properties(std::vector<bool> bounce_back_dir_, 
                                      std::vector<double> bounce_back_delta_);

    /**
     * @brief Set the distribution for the adjacent node alorng direction i.
     * 
     * @param i Direction index
     * @param value 
     */
    inline void set_f_adj(int i, double value) { (*f_adj)[i] = value; }
    
    /**
     * @brief Set the outlet properties object
     * 
     */
    void set_outlet_properties();
    
    /// @}



    /// @name Static lattice variables "D2Q9 model"
    /// @{

    /// Number of dimensions
    static const int dim;
    
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
    static const int dir;

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

    /// Local drag
    // double drag;
    /// Local lift
    // double lift;

    /// Node type (fluid, boundary, solid)
    NodeType node_type;
    /// Directions in which the BCs are to be applied around the node (if it's a boundary node or outlet node)
    std::vector<bool> bounce_back_dir;
    /// Distance between the node and the wall (if it's a boundary node), along the direction of bounce_back_dir
    std::vector<double> bounce_back_delta;
    /// Exit direction for otflow nodes 
    std::vector<bool> outlet_dir;

    /// 2D coordinates of the node
    std::vector<unsigned int> coord;


};




#endif // __NODE_HPP__