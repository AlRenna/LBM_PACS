#ifndef __NODE_HPP__
#define __NODE_HPP__

#include <vector>
class Lattice;

enum class NodeType
{
  fluid,
  boundary,
  solid
};

class Node:
{
  public:
    Node(NodeType node_type_, std::vector<unsigned int> coord_,
        double ux_, double uy_, double rho_);
    ~Node();

    
    

    void init();
    void load_adjacent_velocity_distributions(const Lattice &lattice);
    
    void apply_bc()
    void compute_physical_quantities();

    void collide_stream(const Lattice &lattice);



    // Getters 

    inline double get_f(int i) const { return f[i]; }

    // Setters
    inline void set_boundary_node_dir(std::vector<bool> boundary_node_dir_) { boundary_node_dir = boundary_node_dir_; }




    // Static memebers
    static int dim  = 2;
    static int dir = 9;
    static std::vector<double> weights = {4./9.,
                                         1./9., 1./9., 1./9., 1./9.,
                                         1./36., 1./36., 1./36., 1./36.};
    static std::vector<std::vector<double>> coeff = {{0., 0.},
                                                    {1., 0.}, {-1., 0.}, {0., 1.}, {0., -1.},
                                                    {1., 1.}, {-1., -1.}, {-1., 1.}, {1., -1.}};
    static std::vector<int> bb_indexes = {0, 3, 4, 1, 2, 7, 8, 5, 6};


  private:
    // direction numbering scheme
    // 6 2 5
    // 3 0 1
    // 7 4 8

    // velocitiy distributions
    std::vector<double> f;
    std::vector<double> f_adj;

    // physical quantities
    double ux;
    double uy;
    double rho;

    // node informations
    NodeType node_type;
    std::vector<bool> boundary_node_dir{dir, false}; // in which directions bounce back will be applied (if is a boundary node)
    std::vector<unsigned int> coord;


}




#endif // __NODE_HPP__