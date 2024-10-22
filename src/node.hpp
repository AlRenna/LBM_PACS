/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#ifndef __NODE_HPP__
#define __NODE_HPP__

#include <vector>
class Lattice;

/**
 * @brief Enumeration class for the node types. Used to handle the different Initialization and Boundary conditions.
 * 
 */
enum class NodeType
{
  fluid,
  boundary,
  solid
};

/**
 * @brief Node class. Contains the information of the node in the lattice.
 * 
 */
class Node:
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
    ~Node();

    
    
    /**
     * @brief Initialize the distribution functions to the equilibrium values.
     */
    void init();

    /**
     * @brief Load the adjacent velocity distributions from the adjacent nodes inside the lattice.
     * 
     * @param lattice 
     */
    void load_adjacent_velocity_distributions(const Lattice &lattice);
    
    /**
     * @brief Apply the boundary conditions.
     */
    void apply_bc()

    /**
     * @brief Compute the physical quantities (velocity components and density) from the distribution functions.
     */
    void compute_physical_quantities();

    /**
     * @brief Collision and streaming steps.
     * Compute and relax to the equilibrium distribution. The BGK collision operator is used: 
     * f_i(t+dt) = f_i(t) - 1/tau (f_i(t) - f_i^eq(t))
     * 
     * @param lattice 
     */
    void collide_stream(const Lattice &lattice);



    /// @name Getters
    /// @{ 
    
    /**
     * @brief Get the distribution function at index (direction) i.
     */
    inline double get_f(int i) const { return f[i]; }
    /// @}

    /// @name Setters
    /// @{

    /**
     * @brief Fill the object containg the direction in which the BCs will be applied.
     * 
     * @param boundary_node_dir_ 
     */
    inline void set_boundary_node_dir(std::vector<bool> boundary_node_dir_) { boundary_node_dir = boundary_node_dir_; }
    /// @}



    /// @name Static lattice variables "D2Q9 model"
    /// @{

    /// Number of dimensions
    static int dim  = 2;
    
    /**
     * @brief Number of velocity directions
     * Direction numbering scheme:
     * 6 2 5
     * 3 0 1
     * 7 4 8
     */
    static int dir = 9;

    /// Weights
    static std::vector<double> weights = {4./9.,
                                         1./9., 1./9., 1./9., 1./9.,
                                         1./36., 1./36., 1./36., 1./36.};

    /// Coefficients
    static std::vector<std::vector<double>> coeff = {{0., 0.},
                                                    {1., 0.}, {-1., 0.}, {0., 1.}, {0., -1.},
                                                    {1., 1.}, {-1., -1.}, {-1., 1.}, {1., -1.}};
    
    /// Opposite direction indexes for Bounce Back
    static std::vector<int> bb_indexes = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    /// @}

  private:

    /**
     * @brief Velocity distribution functions.
     * Direction numbering scheme:
     * 6 2 5
     * 3 0 1
     * 7 4 8
     */
    std::vector<double> f;
    /// Velocity distribution functions from adjacent nodes
    std::vector<double> f_adj;

    /// Node velocity components (x)
    double ux;
    /// Node velocity components (y)
    double uy;
    /// Node density
    double rho;

    /// Node type (fluid, boundary, solid)
    NodeType node_type;
    /// Direction in which the BCs are to be applied around the node (if it's a boundary node)
    std::vector<bool> boundary_node_dir{dir, false}; 
    /// 2D coordinates of the node
    std::vector<unsigned int> coord;


}




#endif // __NODE_HPP__