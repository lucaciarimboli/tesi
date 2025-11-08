# Free Boundary Grad–Shafranov Solver

This project implements a solver for the **free boundary Grad–Shafranov problem**, used to compute axisymmetric MHD equilibria in tokamaks. The code is written in Python and uses **Firedrake**, an open-source finite element library.

The current implementation supports both **limiter** and **divertor** configurations and can handle different tokamak geometries. Some extensions are under development, such as imposing physically accurate Neumann boundary conditions via FEM–BEM coupling.

## Project Structure

The solver is organized in an object-oriented fashion. The core class **`GradShafranovSolver`** (in `src/core/solver.py`) stores the problem data and final solution and manages the auxiliary classes which are needed to compute the equilibrium solution.

Auxiliary classes are located in `src/utils`:
- **`Picard`** (`fixed_point.py`): performs fixed-point iterations
- **`Newton`** (`newton.py`): performs Newton iterations
- **`Plasma`** (`plasma.py`): updates the plasma boundary at each iteration
- **`Plot`** (`plot.py`): visualizes the mesh and the computed equilibrium

```
.
|-- src
|   |-- core/solver.py
|   |-- utils
|   |   |-- fixed_point.py
|   |   |-- newton.py
|   |   |-- plasma.py
|   |   |-- boundary_conditions.py
|   |   `-- functions
|   |       |-- mesh_tags.py
|   |       |-- coils.py
|   |       |-- mask.py
|   |       `-- boundary_integrals.py
|   `-- main/main.py
|-- meshes
|-- results
|-- notebooks
`-- BVP
    |-- GradShafranov.ipynb
    `-- meshes
```

The directories `src/`, `meshes/` and `results/` contain all components required to set up the problem, run the solver and inspect the output.

The `notebooks/` folder contains Jupyter notebooks used for testing and prototyping new features. The `BVP/` folder includes an earlier version of the problem where the plasma boundary was fixed. This boundary value problem was solved for multiple plasma shapes and current density profiles, together with a spatial convergence study.

## How to Use

The solver requires the following input data:

- Tokamak geometry and mesh resolution in the vacuum vessel
- Plasma toroidal current density profile \(J_\phi(R, \psi_N)\)
- Its derivative with respect to the normalized flux, needed for Newton iterations
- Currents (in Ampere) in each PF coil

These parameters are defined in the `params` dictionary inside `src/main/main.py`.

By default, the solver starts from a uniform initial flux, performs two Picard iterations to generate an initial guess and then computes the equilibrium via Newton iterations. This behavior, together with iteration parameters and tolerances, can be customized in the same file.

## Tokamak Geometry

Meshes are generated using **Gmsh**. The directory `meshes/` includes an ITER-like geometry (`ITER.geo` and `ITER.msh`).

To use a different geometry, make sure the mesh satisfies these requirements:

1. The limiter is represented as a closed curve embedded in the mesh.
2. Each PF coil is a separate region, with its boundary embedded.
3. Optional passive structures (e.g. vacuum vessel) are also embedded.
4. The external computational boundary is a rectangle; one side must lie on the symmetry axis \(R=0\). The domain should be large enough for homogeneous Dirichlet boundary conditions to be a good approximation.
5. Mesh tags for all regions must be defined in `src/utils/functions/mesh_tags.py`.
6. Coordinates defining limiter, coils and passive structures should be added to `src/utils/plot.py` for correct visualization.

## Further Developments

Two major improvements are planned:

1. **FEM–BEM coupling** to impose physically correct boundary conditions while solving on a reduced computational domain, improving accuracy and reducing degrees of freedom.
2. **Adaptive mesh refinement** inside the vacuum chamber, to improve the accuracy of the X-point and plasma boundary detection while substantially reducing the number of degrees of freedom.

Note. Firedrake does not currently support weak formulations involving double integrals, which appear naturally in FEM–BEM coupling. The class `JN_coupling_BCs` (`src/utils/boundary_conditions.py`) is under development to solve approximately the boundary integral equation appearing in the *Johnson–Nédélec coupling* using a boundary element method, but avoiding direct FEM-BEM coupling.