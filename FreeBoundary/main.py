import sys
sys.path.append(".")

from solver import GradShafranovSolver
from firedrake import Constant

# Problem parameters (adapted from your notebook)
params = {
    "geometry": "build",
    "x0": 0.7,
    "y0": 0.5,
    "R": 0.25,
    "thickness": 0.03,
    "coils": [
        (0.15, 0.25, 0.45, 0.55),  # Coil 1
        (0.25, 0.35, 0.7, 0.8),    # Coil 2
        (0.25, 0.35, 0.2, 0.3),    # Coil 3
    ],
    "domain_size": (1.0, 1.0),
    "mesh_size_min": 0.01,
    "mesh_size_max": 0.05,
    "vacuum_mesh_size": 0.02,
    "limiter_pts": [(0.8, 0.5), (0.6, 0.5)],
    "vessel": (0.7, 0.5, 0.25, 0.03),  # (x0, y0, r, thickness)
    "I": [0.1, 0.1, 0.0],            # Coil currents
    "j_cv": 0.3,                      # Vessel wall current density
    "function_space_family": "P",
    "function_space_degree": 2,
    "max_iterations": 100,
    "tolerance": 1e-5,
    "verbose": True,
    # G function as a lambda (x, psi) -> x**2 + psi
    "G": lambda x, psi: x**2 + psi,
    # Initial guess (can be a Constant or a Firedrake Function)
    "initial_guess": Constant(0.0),
}

if __name__ == "__main__":
    solver = GradShafranovSolver(params)
    solver.display_mesh()
    solver.solve()
    solver.plot_flux()