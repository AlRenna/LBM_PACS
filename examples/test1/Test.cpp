/**
 * @file
 *
 * @author Alessandro Renna <alessandro1.renna@mail.polimi.it>
 * @author Mattia Marzotto <mattia.marzotto@mail.polimi.it>
 */

#include "../../src/lattice.hpp"
#include "../../src/utils_python.cpp"


int main(int argc, char **argv)
{
  // TODO: Modificare script python per far si che i valori di nx e ny vengano presi dal file cpp o da un file di testo
  // TODO: creare un file cpp o python che generi un file .txt contente i valori iniziali e altri parametri necessari all'esecuzione
  // Run the python script to generate the lattice
  run_terminal_python();//{"source env/bin/activate","python3 generate_lattice_RGB.py", "deactivate" });

  int nx = 10;
  int ny = 10;
  // Create the lattice
  Lattice lattice(10, 10, 0.1);

  std::vector<double> ux_in(nx*ny, 0.0);
  std::vector<double> uy_in(nx*ny, 0.0);
  std::vector<double> rho_in(nx*ny, 1.0);

  for (int i = 0; i < nx; ++i)
  {
    ux_in[i] = 0.1;
  }
  // Set the initial and boundary conditions
  lattice.load_ICs_and_BCs(ux_in, uy_in, rho_in, "output.csv");

  // Run the simulation
  lattice.run();

  return 0;
}