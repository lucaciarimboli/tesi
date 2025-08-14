import sys
import time
sys.path.append(".")

from src.core.solver import GradShafranovSolver
from firedrake import Constant
import numpy as np

#--------------------------------------------------#
#                CURRENT PROFILE                   #
#--------------------------------------------------#

# Plasma current density profile G:
def G(R, psi_norm):
    r0 = 6.2
    alpha = 2
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * abs(1 - psi_norm**alpha) ** gamma

# Derivative of G w.r.t. psiN for Newton iterations:
def dGdpsiN(R, psiN):
    r0 = 6.2
    alpha = 2
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * gamma * abs(1 - psiN**alpha) ** (gamma-1) * (-alpha*psiN**(alpha-1))

#--------------------------------------------------#
#                 SET PARAMETERS                   #
#--------------------------------------------------#

# Limiter configurations:
#'I': [-6.705e5, 1.373e4, 2.133e6, 1.432e6, -3.774e5, -6.172e5, -1.885e6, -2.359e6, -2.124e6, -1.836e6, -3.491e6, -2.04e6],
#'I': [-8.208e5, -8.621e4, 2.783e6, 1.703e6, -6.491e5, -7.971e5, -2.026e6, -2.508e6, -2.15e6, -1.874e6, -3.607e6, -2.303e6],

# Divertor configuration:
#'I': [-1.4e6, -9.5e6, -2.04e7, -2.04e7, -1.0e7, 3.6e6, 5.5e6, -2.3e6, -6.5e6, -4.8e6, -7.5e6, 1.73e7], --> non funzion
#'I': [-2.848113e3,-2.205664e+04,-3.022037e4,-3.022037e4,-2.478694e4,1.143284e3,-4.552585e6,3.180596e6,5.678096e6,3.825538e6,1.066498e7,-2.094771e7],   # da paper serino

params = {
    # Tokamak geometry:
    'geometry': "ITER",

    # Currents configuration:
    'I': [-1.4e6, -9.5e6, -2.04e7, -2.04e7, -1.0e7, 3.6e6, 5.5e6, -2.3e6, -6.5e6, -4.8e6, -7.5e6, 1.73e7],
    'j_plasma': G,
    'j_plasma_derivative': dGdpsiN, # needed only for Newton iterations

    # Algorithm parameters:
    'max_iterations': 1000,
    'tolerance': 1e-10,
    'verbose': True,
    'show_plots': True,
    'algorithm': "Newton",

    # Initial guess for the flux function:
    'initial_guess': Constant(0.01) # to avoid infinite normalized error at first iteration
}

#--------------------------------------------------#
#        EXECUTE SOLVER AND PLOT RESULTS           #
#--------------------------------------------------#

if __name__ == "__main__":
    
    solver = GradShafranovSolver(params)
    solver.plot_mesh()
    start_time = time.time()

    #solver.set_algorithm("Picard")
    solver.set_iterations_params(1, params['tolerance'], params['verbose'])

    for i in range(30):
        solver.solve()
        path = f"./results/test2/flux_plot_{i+1}.png"
        solver.plot_flux(path)

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Simulation ended in {minutes} min. and {seconds} sec.")
    #solver.plot_flux()