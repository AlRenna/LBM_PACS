# LBM

## Lattice (.hpp and .cpp)
- Generate and populate the **lattice**;
- run the simulation over each **node**;
- save the output.

## Node (.hpp and .cpp)
- Calculate **collision** and **stream** the information;
- Apply **BounceBack** (Simple and Interpolated);
- Evaluate physical quantities (velocity and density/pressure) and integrals (drag and lift)

## Utils
Contains code to supplement the simulation or other part of the code (poth for **c++** and **Python**)