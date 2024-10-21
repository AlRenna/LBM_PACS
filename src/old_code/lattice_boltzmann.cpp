#include "src/lattice_boltzmann.cpp"

namespace lbm_pacs
{
  lbm::lbm( 
           const double nx_, 
           const double ny_,
            std::shared_ptr<std::vector<double>> rho_,
            std::shared_ptr<std::vector<double>> ux_,
            std::shared_ptr<std::vector<double>> uy_)
    : nx(nx_),
      ny(ny_),
      f0(std::make_shared<std::vector<double>>(nx_*ny_, 0.)),
      f1(std::make_shared<std::vector<double>>(nx_*ny_*8, 0.)),
      f2(std::make_shared<std::vector<double>>(nx_*ny_*8, 0.)),
      rho(rho_),
      ux(ux_),
      uy(uy_)
  {}

  lbm::run()
  {
    init();

    double t = 0.;

    std::make_shared<std::vector<double>> f_tmp;

    for(unsigned int i = 0; i < timesteps; ++i)
    {
      // TODO: sistema struttura
      stream_collide();
      impose_BB_BC();
      f_tmp = f1;
      f1 = f2;
      f2 = f_tmp;

      t = t + dt;
    }
  }
  

  lbm::init()
  {
    /// Initialize the velocity distribution function evaluating the equilibrium distribution

    // TODO: parallelize on gpu
    for(unsigned int y = 0; y < ny; ++y)
    {
        for(unsigned int x = 0; x < nx; ++x)
        {
            double r = rho[scalar_index(x,y)];
            double u  = ux[scalar_index(x,y)];
            double v  = uy[scalar_index(x,y)];
            
            // eval the equilibrium
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]
            
            // temporary variables
            std::vector<double> w_rho = weights*r;
            double omusq = 1.0 - 1.5*(u*u+v*v);
            
            double tu = 3.0*u;
            double tv = 3.0*v;

            f0[scalar_index(x,y)]    = w_rho[0]*(omusq);
            
            for(unsigned int i = 1; i<9; ++i){
              f1[field_index(x,y,i)] = w_rho[i]*(omusq + (coeff[i][0]*tu + coeff[i][1]*tv)*(1.0+0.5*(coeff[i][0]*tu + coeff[i][1]*tv)));
            }
        }

        init_BC();
    }
  }


  lbm::init_BC()
  {
    /// Initialize the map of boundary conditions
    // top, bottom, left, right, 
    // top_left, top_right, bottom_left, bottom_right
    // inlet, outlet, obstacle

    // corners
    BC_map[std::make_pair(0,0)] = "top_left";
    BC_map[std::make_pair(nx-1,0)] = "top_right";
    BC_map[std::make_pair(0,ny-1)] = "bottom_left";
    BC_map[std::make_pair(nx-1,ny-1)] = "bottom_right";

    // sides
    for(unsigned int x = 1; x < nx-1; ++x)
    {
        BC_map[std::make_pair(x,0)] = "top";
        BC_map[std::make_pair(x,ny-1)] = "bottom";
    }
    for(unsigned int y = 1; y < ny-1; ++y)
    {
        BC_map[std::make_pair(0,y)] = "left";
        BC_map[std::make_pair(nx-1,y)] = "right";
    }

    // inlet and outlet

    // obstacle
  }


  lbm::stream_collide()
  {
    /// Streaming and collision
    /// compute the next time step distribution functions along each node

    // TODO: parallelize on gpu
    const double tauinv = 1./tau; // 1/tau
    const double omtauinv = 1.0-tauinv;     // 1 - 1/tau

    std::vector<double> ft{9,0.};

    for(unsigned int y = 0; y < ny; ++y)
    {
        for(unsigned int x = 0; x < nx; ++x)
        {
            // find the adjacent nodes
            unsigned int x_dx = (x+1)%nx;
            unsigned int y_down = (y+1)%ny;
            unsigned int x_sx = (nx+x-1)%nx;
            unsigned int y_up = (ny+y-1)%ny;
            
            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8
            
            ft[0] = f0[scalar_index(x,y)];

            // load populations from adjacent nodes
            ft[1] = f1[field_index(x_sx, y,  1)];
            ft[2] = f1[field_index(x,    y_up,2)];
            ft[3] = f1[field_index(x_dx, y,  3)];
            ft[4] = f1[field_index(x,    y_down,4)];
            ft[5] = f1[field_index(x_sx, y_up,5)];
            ft[6] = f1[field_index(x_dx, y_up,6)];
            ft[7] = f1[field_index(x_dx, y_down,7)];
            ft[8] = f1[field_index(x_sx, y_down,8)];
            
            // compute moments
            double rho_ = 0.0;
            for (unsigned int i = 0; i < 9; ++i)
            {
                rho_ += ft[i];
            }
            double rhoinv = 1.0/rho_;

            double u = 0.;
            double v = 0.;
            for (unsigned int i = 0; i < 9; ++i)
            {
                u += rhoinv*coeff[i][0]*ft[i];
                v += rhoinv*coeff[i][1]*ft[i];
            }
            
            // only write to memory when needed
            // if(save)
            // {
            //     r[scalar_index(x,y)] = rho;
            //     u[scalar_index(x,y)] = ux;
            //     v[scalar_index(x,y)] = uy;
            // }
            
            // Collision and streaming step: now compute and relax to equilibrium 
            // BGK collision operator
            // f_i(t+dt) = f_i(t) - 1/tau (f_i(t) - f_i^eq(t))

            // temporary variables
            std::vector<double> w_tau_rho = weights*rho_*tauinv;
            double omusq = 1.0 - 1.5*(u*u+v*v); // 1-(3/2)u.u
            
            double tu = 3.0*u;
            double tv = 3.0*v;
            
            f0[scalar_index(x,y)]    = omtauinv*ft[0] + w_tau_rho[0]*(omusq);
            
            for(unsigned int i = 1; i<9; ++i){
              f2[field_index(x,y,i)] = omtauinv*ft[i] +
                                      w_tau_rho[i]*(omusq + (coeff[i][0]*tu + coeff[i][1]*tv)*(1.0+0.5*(coeff[i][0]*tu + coeff[i][1]*tv)));
            }
        }
    }
  }

  // TODO: 
  lbm::impose_BB_BC()
  {
    /// Impose the bounce back boundary conditions
    // Read through the map of boundary conditions and apply the bounce back rule
   
    

  }
}