# Free Boundary Grad-Shafranov Solver

This project implements a solver for the free boundary Grad-Shafranov problem, which is commonly used in plasma physics to describe the equilibrium of magnetically confined plasmas in tokamaks.
The solver utilizes a fixed-point iteration method.

## Project Structure

```
freeboundary
├── freeboundary
│   ├── __init__.py
│   ├── grad_shafranov.py
│   ├── gs_varf.py
│   ├── update_plasma_boundary.py
│   ├── initialize_plasma_boundary.py
│   └── coils_contribution.py
├── main.py
└── README.md
```

## Installation

To run this project, you need to have Python installed along with the following dependencies:

- Firedrake
- NumPy
- Matplotlib


## Usage

1. Clone the repository or download the project files.
2. Navigate to the project directory.
3. Run the main script:

```bash
python main.py
```

## Functionality

[...]

## Example

An example of how to set up and solve the Grad-Shafranov problem is provided in the `main.py` file. You can modify the parameters such as coil positions, currents, and vessel geometry to explore different configurations.