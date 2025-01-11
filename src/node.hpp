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
 * @brief Enumeration class for the node types. Used to handle different initialization and boundary conditions.
 * 
 */
enum class NodeType
{
  fluid,
  obstacle,
  solid,
  inlet,
  outlet,
  boundary
};

/**
 * @brief Enumeration class for the ZouHe BCs. Used to handle different positional boundary conditions.
 * 
 */
enum class ZouHeType
{
  none,
  right,
  top,
  left,
  bottom
};

/**
 * @brief Node class. Contains information about the node in the lattice such as type, position, velocity
 * distribution, and physical quantities.
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

    /**
     * @brief Default constructor
     * 
     */
    Node() = default;

    /**
     * @brief Initialize the distribution functions to the equilibrium values.
     */
    void init_equilibrium();

    /**
     * @brief Collision step.
     * 
     * Compute and relax to the equilibrium distribution. The BGK collision operator is used: 
     * 
     * \f$ f_i(t+dt) = f_i(t) - \frac{1}{\tau} (f_i(t) - f_i^{eq}(t)) \f$
     * 
     */
    void collide(const Lattice &lattice);

    /**
     * @brief Streaming step.
     * 
     * Stream the distribution functions to the adjacent nodes (unless BCs are to be applied).
     * 
     */
    void stream(Lattice& lattice);

    /**
     * @brief Apply the Boundary Conditions.
     * 
     */
    void apply_BCs(const Lattice &lattice);

     /**
     * @brief Apply the Interpolated Bounce-Back.
     * 
     * Use the post-collision distribution to compute the IBB:
     * 
     * \f$ f_{bb_i} = 2q_i f_i(x_b) + (1 - 2q) f_i(x_f) \quad \text{if} \quad q_i < \frac{1}{2} \f$
     * 
     * \f$ f_{bb_i} = \frac{1}{2q_i} f_i(x_b) + \left( \frac{2q - 1}{2q} \right) f_{bb_i}(x_b) \quad \text{if} \quad q_i \geq \frac{1}{2} \f$
     *
     */
    void apply_IBB(const Lattice &lattice, unsigned int i);

    /**
     * @brief Apply the Simple Bounce-Back.
     * 
     * Use the post-collision distribution to compute the BB:
     * 
     * \f$ f_{bb_i} = f_i(x_b) \f$
     */
    void apply_BB(const Lattice &lattice, unsigned int i);
    
    /**
     * @brief Apply the ZouHe BCs on outlet nodes by imposing a pressure.
     * 
     */
    void apply_ZouHe(const Lattice &lattice);

    /**
     * @brief Update the distribution functions. Swap the pre-collision with the adjacent distribution pointer for the next iteration.
     * 
     */
    void update_f();

    /**
     * @brief Compute the physical quantities (velocity components and density) from the distribution functions.
     */
    void compute_physical_quantities();

    /**
     * @brief Function to compute physical integrals such as lift and drag.
     * 
     */
    void compute_drag_and_lift(const Lattice &lattice);

    /// @name Getters
    /// @{

    /// Get the horizontal velocity component
    inline double get_ux() const { return ux; }
    /// Get the vertical velocity component
    inline double get_uy() const { return uy; }
    /// Get the density
    inline double get_rho() const { return rho; }
    /// Get the drag value
    inline double get_drag() const { return drag; }
    /// Get the lift value
    inline double get_lift() const { return lift; }
    /// Get the NodeType
    inline NodeType get_node_type() const { return node_type; }
    /// Get the ZouHeType
    inline ZouHeType get_zou_he_type() const { return zou_he_type; }

    /// Get distribution functions pre-collision
    inline const std::vector<double>& get_f_pre() const { return *f_pre; }
    /// Get distribution functions post-collision
    inline const std::vector<double>& get_f_post() const { return *f_post; }
    /// Get adjacent distribution functions 
    inline const std::vector<double>& get_f_adj() const { return *f_adj; }
    /// Get the bounce_back_dir vector
    inline const std::vector<bool>& get_bounce_back_dir() const { return bounce_back_dir; }
    /// Get the bounce_back_delta vector
    inline const std::vector<double>& get_bounce_back_delta() const { return bounce_back_delta; }
    
    /// Get the 2D coordinates of the node
    inline const std::vector<unsigned int>& get_coord() const { return coord; }
    /// @}

    /// @name Setters
    /// @{

    /**
     * @brief Fill the vectors containing the boundary node properties, 
     * such as the direction in which to apply the BCs (IBB) and the 
     * weighted distance between the node and the wall.
     * 
     */
    void set_bounce_back_properties(std::vector<bool> bounce_back_dir_, 
                                      std::vector<double> bounce_back_delta_);

    /**
     * @brief Set the ZouHe type object.
     * 
     */
    void set_zou_he_type(Lattice &lattice);

    /**
     * @brief Set the distribution for the adjacent node along direction i.
     * 
     */
    inline void set_f_adj(int i, double value) { (*f_adj)[i] = value; }
    
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
    
    /// Velocity distribution functions pre-collision
    std::unique_ptr<std::vector<double>> f_pre;
    /// Velocity distribution functions post-collision
    std::unique_ptr<std::vector<double>> f_post;
    /// Velocity distribution functions from adjacent nodes
    std::unique_ptr<std::vector<double>> f_adj;

    /// Node velocity components (x)
    double ux;
    /// Node velocity components (y)
    double uy;
    /// Node density
    double rho;

    /// Local drag
    double drag;
    /// Local lift
    double lift;

    /// Node type (fluid, boundary, solid)
    NodeType node_type;

    /// ZouHeType 
    ZouHeType zou_he_type;

    /// Directions in which the BCs are to be applied around the node (if it's a boundary node or outlet node)
    std::vector<bool> bounce_back_dir;
    /// Distance between the node and the wall (if it's a boundary node), along the direction of bounce_back_dir
    std::vector<double> bounce_back_delta;

    /// 2D coordinates of the node
    std::vector<unsigned int> coord;


};




#endif // __NODE_HPP__