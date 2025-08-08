import sys
sys.path.append(".")

from solver import GradShafranovSolver
from firedrake import *
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

params = {
    # Tokamak geometry:
    "geometry": "ITER",

    # Currents configuration:
    "I": [-8.208e5, -8.621e4, 2.783e6, 1.703e6, -6.491e5, -7.971e5, -2.026e6, -2.508e6, -2.15e6, -1.874e6, -3.607e6, -2.303e6],
    "j_cv": 0.0,               
    "j_plasma": G,
    "j_plasma_derivative": dGdpsiN, # needed only for Newton iterations

    # Algorithm parameters:
    "max_iterations": 1000,
    "tolerance": 1e-5,
    "verbose": True,
    "algorithm": "Picard",

    # Initial guess for the flux function:
    "initial_guess": Constant(0.0)
}

#--------------------------------------------------#
#            TO SHOW GRID INDEPENDENCE             #
#--------------------------------------------------#

if __name__ == "__main__":
    solver = GradShafranovSolver(params)

    psi_max = []
    psi0 = []

    h = [0.4, 0.2, 0.1, 0.05, 0.025]
    # to extract mesh size: m.cell_sizes.dat.data.min()

    for i in range(0,5):
        print('Importing mesh...')
        path = "./meshes/ITER_convergence/ITER" + str(i+1) + ".msh"
        solver.Mesh = Mesh(path, dim = 2, distribution_parameters={"partition": False}, reorder = True)

        print(f'\nSolving with limiter cell size h = {h[i]}...')

        solver.function_spaces()
        solver.BCs = DirichletBC(solver.V, 0.0, solver.tags['boundary'])
        solver.limiter = DirichletBC(solver.V, 0.0, solver.tags['limiter']).nodes
        
        #solver.set_iterations_params(max_iterations=50, tolerance=params["tolerance"], verbose=True)
        #solver.set_algorithm("Picard")
        #solver.solve()

        #solver.set_iterations_params(max_iterations=1000, tolerance=params["tolerance"], verbose=True)
        #solver.set_algorithm("Marder-Weitzner")
        #solver.set_initial_guess(initial_guess=solver.psi)
        
        solver.solve()

        psi_max.append(solver.psi_max)
        psi0.append(solver.psi0)

    psi_max = np.array(psi_max) / max(psi_max)
    psi0 = np.array(psi0) / max(psi0)
    #h_1 = np.array([1.0 / hi for hi in h]) # h.^(-1)
    h = np.array(h)
    
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(h, psi_max, marker='o', label='psi_max')
    plt.plot(h, psi0, marker='s', label='psi0')
    plt.xlabel('h')
    plt.ylabel('Value')
    plt.title('Picard Convergence:')
    plt.legend()
    plt.grid(True)
    plt.savefig("./results/grid_independence.png")
