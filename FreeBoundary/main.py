import sys
import time
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

lim_line = circle_points(x0, y0, radius - 2*thick, n=100)

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
#                CURRENT PROFILE                   #
#--------------------------------------------------#

def G(R, psi_norm):
    r0 = 6.2
    alpha = 2.0
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * (1 - psi_norm**alpha) ** gamma

#--------------------------------------------------#
#                 SET PARAMETERS                   #
#--------------------------------------------------#

params = {

    "geometry": "ITER",
    #"boundary": default,
    #"vessel_outer_wall": vessel_outer,
    #"vessel_inner_wall": vessel_inner,
    #"coils": [coil1, coil2],
    #"limiter_pts": None,
    #"limiter_line": lim_line,
    #"mesh_size_min": 0.005,
    #"mesh_size_max": 0.01,
    #"limiter_mesh_size": 0.001,
    #"dist_from_limiter": 0.1,
    #"coils_mesh_size": 0.01,
    #"dist_from_coils": 0.1,
    "I": [-6.705e5, 1.373e4, 2.133e6, 1.432e6, -3.774e5, -6.172e5, -1.885e6, -2.359e6, -2.124e6, -1.836e6, -3.491e6, -2.04e6],
    "j_cv": 0,                # Vessel wall current density
    "function_space_family": "P",
    "function_space_degree": 2,
    "max_iterations": 1000,
    "tolerance": 1e-4,
    "verbose": True,
    "G": G,
    # Initial guess (can be a Constant or a Firedrake Function)
    "initial_guess": Constant(0.0),
    "algorithm": "Picard",
    #"algorithm": "Marder-Weitzner",
    #"alpha": 0.4,  # Relaxation Parameter
    #"algorithm": "Newton"
}


#--------------------------------------------------#
#        EXECUTE SOLVER AND PLOT RESULTS           #
#--------------------------------------------------#


if __name__ == "__main__":
    
    start_time = time.time()

    params["tolerance"] = 1e-2
    params["algorithm"] = "Newton"
    solver = GradShafranovSolver(params)
    #solver.display_mesh()
    solver.solve()
    #solver.plot_flux()

    solver.set_algorithm("Marder-Weitzner")
    solver.solve()

    solver.set_algorithm("Picard")
    solver.solve()


    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Simulation ended in {minutes} min. and {seconds} sec.")


#--------------------------------------------------#
#            TO SHOW GRID INDEPENDENCE             #
#--------------------------------------------------#
'''
if __name__ == "__main__":
    solver = GradShafranovSolver(params)
    solver.display_mesh()
    solver.convergence_test()
    

    psi_max = []
    psi0 = []

    for i in range(6):
        params["mesh_size_min"] = 0.01 / (2 ** i)
        params["mesh_size_max"] = 0.05 / (2 ** i)
        params["limiter_mesh_size"] = 0.1 / (2 ** i)
        params["dist_from_limiter"] = 0.1 / (2 ** i)

        print(f"Running solver with mesh size factor: {2**i}")
        solver = GradShafranovSolver(params)
        solver.display_mesh()
        solver.solve()
        solver.plot_flux()

        psi_max.append(solver.psi.vector().max())
        psi0.append(solver.psi0)
    
    import matplotlib.pyplot as plt

    refinement_levels = list(range(6))
    plt.figure()
    plt.plot(refinement_levels, psi_max, marker='o', label='psi_max')
    plt.plot(refinement_levels, psi0, marker='s', label='psi0')
    plt.xlabel('Refinement level (i)')
    plt.ylabel('Value')
    plt.title('Grid Independence Study')
    plt.legend()
    plt.grid(True)
    plt.savefig("./results/grid_independence.png")
'''

#--------------------------------------------------#
#     SOME PICARD ITERATIONS FOR INITIAL GUESS     #
#--------------------------------------------------#
'''
if __name__ == "__main__":

    # Perform some Picard iterations:
    params["max_iterations"] = 50
    #params["algorithm"] = "Picard"
    params["algorithm"] = "Marder-Weitzner"
    params["alpha"] = 0.4 # relaxation parameter

    solver = GradShafranovSolver(params)
    solver.solve()

    # Use the result as initial guess for Newton method:
    solver.set_iterations_params(max_iterations=1000, tolerance=params["tolerance"], verbose=True)
    solver.set_algorithm("Newton")
    solver.set_initial_guess(initial_guess=solver.psi)

    # Solve using Newton method starting from a closer psi_initial:
    solver.solve()
    solver.plot_flux()
'''