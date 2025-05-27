# Free Boundary Grad-Shafranov Solver

This project implements a solver for the Grad-Shafranov equation, which is commonly used in plasma physics to describe the equilibrium of magnetically confined plasmas. The solver utilizes a fixed-point iteration method to compute the plasma boundary and magnetic flux function.

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
- SymPy

You can install the required packages using pip:

```bash
pip install firedrake numpy matplotlib sympy
```

## Usage

1. Clone the repository or download the project files.
2. Navigate to the project directory.
3. Run the main script:

```bash
python main.py
```

## Functionality

- **GradShafranov**: Solves the Grad-Shafranov equation using the provided parameters, including coil currents, vessel geometry, and initial conditions.
- **Variational Forms**: Defines the variational forms used in the Grad-Shafranov problem.
- **Plasma Boundary Update**: Updates the plasma boundary based on the computed flux function.
- **Coils Contribution**: Computes the contribution of the coils to the variational problem.

## Example

An example of how to set up and solve the Grad-Shafranov problem is provided in the `main.py` file. You can modify the parameters such as coil positions, currents, and vessel geometry to explore different configurations.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.