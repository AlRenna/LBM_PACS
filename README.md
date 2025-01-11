# LBM_PACS

## Index
- [Setup](#setup)
  - [Required packages](#required-packages)
  - [Optional packages](#optional-packages)
  - [Create the environment](#create-the-environment)
  - [Create the build folder](#create-the-build-folder)
  - [Compile and execute an example](#compile-and-execute-an-example)
  - [Generate Code Documentation (Doxygen)](#generate-code-documentation-doxygen)
- [The LBM implementation](#the-lbm-implementation)
  - [Preprocessing](#preprocessing)
  - [Running the simulation](#running-the-simulation)
  - [Postprocessing](#postprocessing)
- [How to build an example](#how-to-build-an-example)
- [References](#references)

## SETUP
### Required packages
- OpenMp
- ffmpeg 

### Optional packages
- CUDA Toolkit (12.6)
- doxygen
- graphviz (for documentation)

### Create the environment
Once you run the following commands a python environment with all of its dependencies will be installed. Remember to always have your environment activated.
```bash
python3 -m venv env
source env/bin/activate
pip install -r py_requirements.txt
```

### Create the build folder 
This folder will contain all of the generated content. Executables will be compiled here and all of the simulation outputs will be stored in this directory. 
```bash
mkdir build
cd build
cmake ..
make
```

### Compile and execute an example
Once you run the `make` command you will find each executable inside the corresponding folder in examples.
```bash
cd build/examples/chose_your_test
./test_executable
```

To run the test with different properties modify the `params.json` file to control the simulation parameters, such as the image of the fluid dynamics environment, grid size, boundary conditions, and  viscosity.

```bash
./test_executable -gpu
```

If you add the `-gpu` argument, the code will run using a GPU parallelization (If CUDA Toolkit is not installed the code will ignore this and run using OpenMP parallelization).

### Generate Code Documentation (Doxygen)
The docs folder contains the documentation of the repository. 
```bash
doxygen Doxyfile
```

## The LBM implementation
The goal of this code is to simulate fluid dynamics using the Lattice Boltzmann Method (LBM). This computational approach is particularly useful for modeling fluid systems with a relatively simple numeric scheme the can be extended for parallel computation. 

### Preprocessing

Once the user runs the executable, the `lattice_generation_RGB.py` file preprocesses the data encoded in the image specified in `params.json`. The script outputs a `.csv` file with the list of nodes and its features. 
In order for the python script to correctly classify the nodes of the lattice there is a specific color coding to assign the correct type according to the pixel color.

![Example Image](/images/channel_obs.png)

**NOTE:**  The values of `nx` and `ny` must be large enough to be able to select the contour of the image.

These are the corresponding RGB encodings:
- **Fluid:** black (0,0,0)
- **Solid:** blue (0,0,255)
- **Inlet:** green (0,255,0)
- **Outlet:** red (255,0,0)
- **Obstacle:** white (255,255,255)


### Running the simulation

After the preprocessing step the executable will create a `Lattice` object and populate it using the information in `lattice.csv` and `params.json`. The simulation starts calling the `run` function.

 For each timestep we loop the following operations: 
- **Collision step:** collide using the BGK operator
- **Streaming step:** stream the post collision distributions into the adjacent nodes. 
- **BCs application:** apply the chosen Bounce-Back method depending on the node type
- **Compute drag and lift:** only applied if an obstacle is present
- **Update distributions:** swap pointers for the pre-collision and adjacent distribution for the ext timestep 
- **Compute physical quantities:** calculates velocity and pressure

**Note:** The collision and streaming step must be done for all nodes before moving forward. 


### Postprocessing
At the end, the `animation.py` script will elaborate the output files of the simulation and create plots and animations. 

## How to build an example

The solver is built to be easily used, the user can add new examples to solve different problems.

Create a `new_test_example` folder inside `examples` and  `new_test_example.cpp`, `CMakeLists.txt` and `params.json` inside it. These files should have the following structures:

### CMakeLists.txt

```bash
cmake_minimum_required(VERSION 3.10)

# Collect all source files in the directory
file(GLOB_RECURSE NEW_TEST_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Create an executable from the source files
add_executable(new_test_example ${NEW_TEST_EXAMPLE_SOURCES})

# Link the executable with the main library
target_link_libraries(new_test_example PRIVATE LatticeLib)

# Copy params.json to the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/params.json DESTINATION ${CMAKE_BINARY_DIR}/examples/new_test_example)
```

### params.json

```bash
"image_path": "../../images/new_test_example_image.png",
"lattice": {
    "nx": 300,
    "ny": 300
},
"time": {
    "T_final": 100.0,
    "save_time": 1.0
},
"test_info": {
    "Length": 10.0,
    "velocity": 1.0,
    "rho": 1.0,
    "nu": 0.01
},
"generated_variables":{
    
} 
```

nx and ny are the lattice size. T_final and save_time represent the total simulation time and the time after which data is saved.
The test_info section contains the system’s characteristic length and speed, density and viscosity.
The section genereted_variables will be filled by the code during execution and will contain the variables converted into.

The user can add more section if a more precise set up is required. In our tests, we have used velocity in test_info as the constant inlet velocity.

### new_test_example.cpp
```bash

#include "src/lattice.hpp"
#include "src/utils/utils_python.hpp"

#include <iostream> 
#include <stdexcept>
#include <fstream>
#include <nlohmann/json.hpp>


int main(int argc, char **argv)
{
    
    // Generate the lattice information and variables
    std::string script_path = "../../src/python_scripts/lattice_generation_RGB.py";
    std::cout << "Running python script: lattice_generation_RGB.py" << std::endl;
    run_terminal_python(script_path);

    // Read params.json
    std::ifstream params_file("params.json");
    nlohmann::json params;
    params_file >> params;
    params_file.close();
    ...

    // Create the lattice
    std::cout << "Creating lattice" << std::endl;
    Lattice lattice;
    const int nx = lattice.get_nx();
    const int ny = lattice.get_ny();
    std::vector<double> ux_in(nx*ny, ...);
    std::vector<double> uy_in(nx*ny, ...);
    std::vector<double> rho_in(nx*ny, ...);

    // Modify the initial conditions according to the problem
    ...

    // Set the initial and boundary conditions
    lattice.initialize(ux_in, uy_in, rho_in, "lattice.csv");

    // Run the simulation
    lattice.run(argc, argv);

    // Generate the animation for velocity and density, and plot for drag and lift
    script_path = "../../src/python_scripts/animation.py";
    std::cout << "Running python script:  animation.py"  << std::endl;
    run_terminal_python(script_path);

    return 0;
}
```

The user should initialize the vectors for the velocity and density according to the new problem. It may use one of the functions inside `utils.cpp` or implement new ones. 

**NOTE:** 
The velocity value will not be used as a constant value during a simulation, but it will be multiplied by the following function:

$$f(t) = \frac{1}{1 + \exp(-25 \cdot (t - 0.2))},$$

for example in the case of an inlet we would have:

$$v_{inlet}(t) = ConstantVelocity * f(t).$$


## References

1. Krüger, T., Kusumaatmaja, H., Kuzmin, A., Shardt, O., Silva, G., & Viggen, E. M. (2017). The Lattice Boltzmann Method: Principles and Practice. Springer.

2. Inamuro, T., Yoshino, M., & Suzuki, K. (2021). An Introduction to the Lattice Boltzmann Method. World Scientific. doi:10.1142/12375 

3. Pacheco, P. S., & Malensek, M. (2020). An Introduction to Parallel Programming (2nd ed.). Morgan Kaufmann.

4. Izquierdo, S., & Fueyo, N. (2008). Characteristic nonreflecting boundary conditions for open boundaries in lattice Boltzmann methods. Physical Review E, 78(4), 046707. doi:10.1103/PhysRevE.78.046707

5. Inamuro, T., Yoshino, M., & Ogino, F. (1995). A non-slip boundary condition for lattice Boltzmann simulations. Physics of Fluids.
