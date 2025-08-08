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
    alpha = 2.0
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * abs(1 - psi_norm**alpha) ** gamma

# Derivative of G w.r.t. psiN for Newton iterations:
def dGdpsiN(R, psiN):
    r0 = 6.2
    alpha = 2.0
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
#'I': [-1.4e6, -9.5e6, -2.04e7, -2.04e7, -1.0e7, 3.6e6, 5.5e6, -2.3e6, -6.5e6, -4.8e6, -7.5e6, 1.73e7],

params = {
    # Tokamak geometry:
    'geometry': "ITER",

    # Currents configuration:
    'I': [-1.4e6, -9.5e6, -2.04e7, -2.04e7, -1.0e7, 3.6e6, 5.5e6, -2.3e6, -6.5e6, -4.8e6, -7.5e6, 1.73e7],
    'j_cv': 0.0,               
    'j_plasma': G,
    'j_plasma_derivative': dGdpsiN, # needed only for Newton iterations

    # Algorithm parameters:
    'max_iterations': 1000,
    'tolerance': 1e-5,
    'verbose': True,
    'show_plots': False,
    'algorithm': "Newton",

    # Initial guess for the flux function:
    'initial_guess': Constant(0.0)
}

#--------------------------------------------------#
#        EXECUTE SOLVER AND PLOT RESULTS           #
#--------------------------------------------------#

if __name__ == "__main__":
    
    start_time = time.time()
    solver = GradShafranovSolver(params)

    #solver.set_algorithm("Picard")
    #solver.set_iterations_params(5, params['tolerance'], params['verbose'])
    solver.solve()

    #solver.set_algorithm("Newton")
    #solver.set_iterations_params(params['max_iterations'], params['tolerance'], params['verbose'])
    #solver.solve()

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Simulation ended in {minutes} min. and {seconds} sec.")