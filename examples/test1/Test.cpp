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
    
    std::cout << "Running python script" << std::endl;
    run_terminal_python(script_path);

    
    // Update nx and ny with the new values returned by the Python script
   

    // Create the lattice
    std::cout << "Creating lattice" << std::endl;
    Lattice lattice;
    const int nx = lattice.get_nx();
    const int ny = lattice.get_ny();

    std::vector<double> ux_in(nx*ny, 0.0);
    std::vector<double> uy_in(nx*ny, 0.0);
    std::vector<double> rho_in(nx*ny, 1.0);

    for (int i = 0; i < nx; ++i) {
        ux_in[i] = 0.1;
    }
    // Set the initial and boundary conditions
    lattice.load_ICs_and_BCs(ux_in, uy_in, rho_in, "output.csv");

    // Run the simulation
    lattice.run();

    return 0;
}