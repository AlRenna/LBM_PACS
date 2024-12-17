/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

// #include <Python.h>
#include "src/lattice.hpp"
#include "src/utils/utils_python.hpp"

#include <iostream> 
#include <stdexcept>
#include <fstream>
#include <nlohmann/json.hpp>

int main(int argc, char **argv)
{
    // Read parameters from JSON file
    std::ifstream params_file("params.json");
    nlohmann::json params;
    params_file >> params;
    params_file.close();
    double u_in = params["test_info"]["u_in"];    
    double rho_fluid = params["test_info"]["rho"];

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
    std::vector<double> rho_in(nx*ny, rho_fluid);

    // Set the inlet velocity on left wall
    ux_in = uniform_left_inlet(u_in, nx, ny, "lattice.csv");

    // Set the initial and boundary conditions
    lattice.initialize(ux_in, uy_in, rho_in, "lattice.csv");

    // Run the simulation
    lattice.run(argc, argv);

    // Path to the python script to postprocess the results and create the animation
    script_path = "../../src/python_scripts/animation.py";
    std::cout << "Running python script:  animation.py"  << std::endl;
    run_terminal_python(script_path);

    return 0;
}