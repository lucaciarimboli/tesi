import sys
sys.path.append(".")

from solver import GradShafranovSolver
from firedrake import Constant
import numpy as np


#--------------------------------------------------#
#        DEFINE SPLINES FOR CUSTOMIZE TOKAMAK      #
#--------------------------------------------------#

# Outer boundary (1.0, 1.0) square domain by default

# Inner and outer vessel walls as circles:
def circle_points(x0, y0, r, n=100):
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return [(x0 + r * np.cos(t), y0 + r * np.sin(t)) for t in theta]

x0 = 0.7
y0 = 0.5
radius = 0.25
thick = 0.03

vessel_outer = circle_points(x0, y0, radius + thick, n=100)
vessel_inner = circle_points(x0, y0, radius, n=100)

# Two squared coils:
def rectangle_points(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]

coils_adapted_format = [
    [0.15, 0.25, 0.45, 0.55],
    [0.25, 0.35, 0.7, 0.8],
]

coil1 = rectangle_points(coils_adapted_format[0][0], coils_adapted_format[0][1],
                         coils_adapted_format[0][2], coils_adapted_format[0][3])
coil2 = rectangle_points(coils_adapted_format[1][0], coils_adapted_format[1][1],
                         coils_adapted_format[1][2], coils_adapted_format[1][3])

#--------------------------------------------------#
#                 SET PARAMETERS                   #
#--------------------------------------------------#

params = {

    #  Da rimuovere dopo aver fixato il problema dei subdomains  #
    "x0": x0,
    "y0": y0,
    "R": radius,
    "thickness": thick,
    "vessel": (x0, y0, radius, thick),
    "coils_adapted_format": coils_adapted_format,
    # ---------------------------------------------------------- #

    "geometry": "build",
    #"boundary": default,
    "vessel_outer_wall": vessel_outer,
    "vessel_inner_wall": vessel_inner,
    "coils": [coil1, coil2],
    "limiter_pts": [(0.8, 0.5), (0.6, 0.5)],
    # "limiter_line": None,
    "mesh_size_min": 0.005,
    "mesh_size_max": 0.01,
    "limiter_mesh_size": 0.001,
    "dist_from_limiter": 0.1,
    #"coils_mesh_size": 0.01,
    #"dist_from_coils": 0.1,
    "I": [0.1, 0.1],        # Coil currents
    "j_cv": 3.0,                # Vessel wall current density
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


#--------------------------------------------------#
#        EXECUTE SOLVER AND PLOT RESULTS           #
#--------------------------------------------------#

if __name__ == "__main__":
    solver = GradShafranovSolver(params)
    solver.display_mesh()
    solver.solve()
    solver.plot_flux()