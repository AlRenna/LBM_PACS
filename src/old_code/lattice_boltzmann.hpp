#ifndef __LB_METHOD_HPP_
#define __LB_METHOD_HPP_

#include "src/lattice.hpp"

#include <string>
#include <vector>


namespace lbm_pacs
{

class lbm
  {
  public:
    lbm(const double nx_, 
        const double ny_,
        std::shared_ptr<std::vector<double>> rho_, // TODO: allocare lo spazio?
        std::shared_ptr<std::vector<double>> ux_,
        std::shared_ptr<std::vector<double>> uy_);

    void run();

  private:
    
    void init();
    void init_BC(); // determina quali nodi sono sui bordi e che tipo di BC applicare

    void stream_collide();
    void impose_BC();


    inline size_t scalar_index(unsigned int x, unsigned int y)
    {
      return nx*y+x;
    }

    inline size_t field_index(unsigned int x, unsigned int y, unsigned int component)
    {
      return (nx*y+x)*8+(component-1); // 8 = num of directions - 1 
    }
    

    // DATA

    const double nx;
    const double ny;

    const std::vector<std::vector<double>> coeff = {{ 0., 0. },
                                              { 1., 0. }, { -1., 0. },
                                              { 0., 1. }, { 0., -1. },
                                              { 1., 1. }, { -1., -1. },
                                              { -1., 1. }, { 1., -1. }};

    const std::vector<double> weights = {4./9.,
                                  1./9., 1./9., 1./9., 1./9.,
                                  1./36., 1./36., 1./36., 1./36.};

    std::shared_ptr<std::vector<double>> f0;
    std::shared_ptr<std::vector<double>> f1;
    std::shared_ptr<std::vector<double>> f2;
    std::shared_ptr<std::vector<double>> rho:
    std::shared_ptr<std::vector<double>> ux;
    std::shared_ptr<std::vector<double>> uy; 
    
    unsigned double dt = 0.1;
    unsigned int timesteps = 100;

    // index for the bounce back BCs (index along the opposite direction)
    const std::vector<unsigned int> bb_index = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    std::map<std::pair<unsigned int, unsigned int>, std::string> BC_map;

    

    // TODO: get value from input
    // Viscosity
    const double nu = 1.0/6.0; 
    // Relaxation time
    const double tau = 3.0*nu+0.5;

  }
}


#endif // __LB_METHOD_HPP_