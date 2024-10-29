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

  // Create the lattice
  Lattice lattice(100, 100, 0.1);

  // Set the initial and boundary conditions
  // lattice.set_ICs_&_BCs({0.1, 0.1}, {0.1, 0.1}, {0.1, 0.1}, {NodeType::FLUID, NodeType::FLUID}, {{true, true}, {true, true}});

  // Run the simulation
  // lattice.run();

  return 0;
}