#ifndef __LATTICE_HPP_
#define __LATTICE_HPP_

#include <vector>

namespace lbm_pacs
{

class lattice
  {
  public:
    lattice(/* args */);
    ~lattice();


  private:
    // 
    const double x_min;
    const double x_max;
    const double y_min;
    const double y_max;
    
    const double nx;
    const double ny;
    const std::vector<double> Coeff;
  };
  
  
}


#endif // __LATTICE_HPP_