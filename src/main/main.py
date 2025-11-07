import sys
import time
sys.path.append(".")

from src.core.solver import GradShafranovSolver
from firedrake import *
from src.utils.functions.CompassU import linear_interpolation, CompassU_J

#--------------------------------------------------#
#                 PLASMA CURRENT                   #
#--------------------------------------------------#
'''
# Plasma current density profile j_phi:
def jphi(R, psi_norm):
    r0 = 6.2
    alpha = 2
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * abs(1 - psi_norm**alpha) ** gamma

# Derivative of j_phi w.r.t. psiN for Newton iterations:
def jphi_derivative(R, psiN):
    r0 = 6.2
    alpha = 2
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * gamma * abs(1 - psiN**alpha) ** (gamma-1) * (-alpha*psiN**(alpha-1))

'''
# FOR COMPASS UPGRADE GEOMETRY:
mu0 = 4e-7 * pi
eq_idx = 0
psiN_data, pprime_data, FFprime_data = CompassU_J(eq_idx)

def pprime(psiN):
    return linear_interpolation(psiN,psiN_data,pprime_data)
def FFprime(psiN):
    return linear_interpolation(psiN,psiN_data,FFprime_data)

def jphi(R, psiN):
    return R * pprime(psiN) + 1/(mu0*R) * FFprime(psiN)
def jphi_derivative(R, psiN):
    return diff(jphi(R, psiN), psiN)


#--------------------------------------------------#
#               PF COILS CURRENTS                  #
#--------------------------------------------------#
'''
# Limiter configurations:
'I': [-6.705e5, 1.373e4, 2.133e6, 1.432e6, -3.774e5, -6.172e5, -1.885e6, -2.359e6, -2.124e6, -1.836e6, -3.491e6, -2.04e6],
'I': [-8.208e5, -8.621e4, 2.783e6, 1.703e6, -6.491e5, -7.971e5, -2.026e6, -2.508e6, -2.15e6, -1.874e6, -3.607e6, -2.303e6],

# Divertor configuration:
'I': [-1.4e6, -9.5e6, -2.04e7, -2.04e7, -1.0e7, 3.6e6, 5.5e6, -2.3e6, -6.5e6, -4.8e6, -7.5e6, 1.73e7],
'''

#--------------------------------------------------#
#                  SET PARAMETERS                  #
#--------------------------------------------------#
params = {
    # Tokamak geometry:
    'geometry': "CompassU", # "ITER" or "CompassU"
    'h': 0.05, # mesh size in the vacuum region

    # Currents configuration:
    'I': [
        -654043.86, -449965.2136, -449965.2136, -63607.206, 459332.916,
        578029.5081, 578029.5081, 164952.4866, 164952.4866, 222340.6416,
        50441.787, -601765.3804, -302472.6912, -487972.6832, -485041.4354
    ],
    'j_plasma': jphi,
    'j_plasma_derivative': jphi_derivative, # needed only for Newton iterations

    # Algorithm parameters:
    'max_iterations': 100,
    'tolerance': 1e-10,
    'verbose': True,
    'show_plots': True,
    'algorithm': "Newton",  # "Picard" or "Newton"

    # Initial guess for the flux function:
    'initial_guess': Constant(0.01) # to avoid infinite normalized error at first iteration
}


#--------------------------------------------------#
#        EXECUTE SOLVER AND PLOT RESULTS           #
#--------------------------------------------------#
if __name__ == "__main__":
    
    # INITIALIZE SOLVER:
    solver = GradShafranovSolver(params)
    #solver.plot_mesh()

    # SET INITIAL CONDITION WITH 2 PICARD ITERATIONS:
    solver.set_algorithm("Picard")
    solver.set_iterations_params(3, params['tolerance'], params['verbose'])
    solver.solve()
    solver.set_algorithm("Newton")
    solver.set_iterations_params(params['max_iterations'], params['tolerance'], params['verbose'])
    #path = f"./results/initial_guess.png"
    #solver.plot_flux(path)
    
    # SOLVE THE PROBLEM:
    start_time = time.time()
    solver.solve()

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Simulation ended in {minutes} min. and {seconds} sec.")
    solver.plot_flux()