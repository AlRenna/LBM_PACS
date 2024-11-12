/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

// #include <Python.h>
#include "../../src/lattice.hpp"
#include "../../src/utils_python.cpp"

#include <iostream> 
#include <stdexcept>

int main(int argc, char **argv)
{
    // Path to the python script to preprocess the image and classify the nodes
    std::string script_path = "../../src/python_scripts/lattice_generation_RGB.py";
    std::cout << "Running python script: lattice_generation_RGB.py" << std::endl;
    run_terminal_python(script_path);
   
    // Create the lattice
    std::cout << "Creating lattice" << std::endl;
    Lattice lattice;
    const int nx = lattice.get_nx();
    const int ny = lattice.get_ny();
    std::vector<double> ux_in(nx*ny, 0.0);
    std::vector<double> uy_in(nx*ny, 0.0);
    std::vector<double> rho_in(nx*ny, 1.0);

    // Set the lid driven boundary condition on the upper wall
    ux_in = lid_driven(10, nx, ny, "output.csv");

    // Set the initial and boundary conditions
    lattice.load_ICs_and_BCs(ux_in, uy_in, rho_in, "output.csv");

    // Run the simulation
    lattice.run();

    // Path to the python script to postprocess the results and create the animation
    script_path = "../../src/python_scripts/animation.py";
    std::cout << "Running python script:  animation.py"  << std::endl;
    run_terminal_python(script_path);

    return 0;
}