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
    lbm(/* args */);

    void run();

    void init();
    void collide();
    void stream();
    void update();

  private:

  std::shared_ptr<lattice> lattice_ptr;

  std::vector<double> f0;
  std::vector<double> f1;
  std::vector<double> f2;
  std::vector<double> rho:
  std::vector<double> ux;
  std::vector<double> uy; 
  
  }
}


#endif // __LB_METHOD_HPP_