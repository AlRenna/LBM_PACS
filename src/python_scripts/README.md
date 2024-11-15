## lattice_generation_RGB.py
- Read the given image, nx and ny
- Create a **lattice** of **nodes** over the image, adpating the number of point (nx and ny) if their are not sufficient to evenly cover the entire surface
- Identify the **Type** of the node by check the color:
  - **Red: Wall**
  - **Green: Inlet**
  - **Blue: Outlet**
  - **Black: Fluid**
- Also identify **Boundary** Nodes (Fluid nodes near **Wall/solid nodes**)
- Return a **.csv** file cointaing the information of the lattice:
  - Nodes coordinate
  - Nodes type
  - Nodes distances to wall nodes (for **Interpolated BounceBack**)


## animation.py
Given the output file of a simulation it generates animations for the velocity field and the density.