# LBM_PACS

## SETUP
### Required packages
- OpenMp
- ffmpeg 
- 
### Optional packages
- CUDA Toolkit (12.6)
- doxygen

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


## References

1. Krüger, T., Kusumaatmaja, H., Kuzmin, A., Shardt, O., Silva, G., & Viggen, E. M. (2017). The Lattice Boltzmann Method: Principles and Practice. Springer.

2. Inamuro, T., Yoshino, M., & Suzuki, K. (2021). An Introduction to the Lattice Boltzmann Method. World Scientific. doi:10.1142/12375 

3. Pacheco, P. S., & Malensek, M. (2020). An Introduction to Parallel Programming (2nd ed.). Morgan Kaufmann.

4. Izquierdo, S., & Fueyo, N. (2008). Characteristic nonreflecting boundary conditions for open boundaries in lattice Boltzmann methods. Physical Review E, 78(4), 046707. doi:10.1103/PhysRevE.78.046707

5. Inamuro, T., Yoshino, M., & Ogino, F. (1995). A non-slip boundary condition for lattice Boltzmann simulations. Physics of Fluids.
