from ufl import grad
from firedrake import *

from src.utils.plasma import Plasma
from src.utils.boundary_conditions import JN_coupling_BCs
from src.utils.fixed_point import Picard
from src.utils.newton import Newton
from src.utils.plot import Plot
from src.utils.functions.mesh_tags import get_tags
from src.utils.functions.coils import compute_j_coils

import os
import matplotlib.pyplot as plt
import numpy as np


class GradShafranovSolver:

    def __init__(self, params):
        """
        Constructor for the GradShafranovSolver class.
        This class is designed to solve the Grad-Shafranov equation using a finite element method.
        It initializes the solver with the given parameters and prepares the mesh for the problem.
        
        @param params: Dictionary containing the parameters for the solver.
        """
        # Identify mesh geometry
        geometry = params.get('geometry', "ITER")

        # Initialize the mesh and function spaces
        self.import_mesh(geometry)
        self.V = FunctionSpace(self.m, "P", 1)

        # Identify the tags for the mesh regions and elements:
        self.tags = get_tags(geometry)

        # Initialize flux function:
        self.psi = Function(self.V).interpolate(params.get('initial_guess', Constant(0.0)))

        # Initialize plasma masks and values
        self.plasma = Plasma(self.V,params['h'],self.tags['limiter'],self.tags['vacuum'])

        # Compute the coils current density:
        j_coils = compute_j_coils(self.m, self.tags['coils'], params["I"])

        # Save input data for the problem in a structure:
        self.data = {
            'G': params['j_plasma'],
            'j_coils': j_coils,
            'j_cv': params.get('j_cv', 0.0),
            'dGdpsiN': params.get('j_plasma_derivative',None)
        }

        # Set boundary conditions
        #self.dirichlet = DirichletBC(self.V, 0.0, self.tags['dirichlet_boundary'])
        #self.neumann = JN_coupling_BCs(self.V, self.tags['neumann_boundary'])
        self.dirichlet = DirichletBC(self.V, 0.0, self.tags['boundary']) # homogeneous

        # Set algorithm parameters:
        self.set_algorithm(params.get('algorithm', "Picard"))
        self.set_iterations_params(params.get('max_iterations',100), params.get('tolerance',1e-5), params.get('verbose',False))

        # Define object for plotting the results:
        if params.get('show_plots',False):
            self.plots_flag = True
            self.plot = Plot(self.m, geometry)
        else:
            self.plots_flag = False

    # SETTERS:
    def set_algorithm(self, algorithm):
        """
        Set the algorithm for solving the Grad-Shafranov equation.

        @param algorithm: The algorithm to use (e.g., "Picard", "Newton").
        """
        if algorithm not in ["Picard", "Newton"]:
            self.algorithm = "Picard"  # Default to Picard if invalid algorithm is provided
            print(f"Invalid algorithm '{algorithm}' provided. Defaulting to Picard (fixed point iterations).")
        else:
            self.algorithm = algorithm

    def set_iterations_params(self, max_iterations, tolerance, verbose):
        """
        Set the parameters for the solver iterations.
        
        @param max_iterations: Maximum number of iterations for the solver.
        @param tolerance: Tolerance for convergence.
        @param verbose: If True, print iteration details.
        """
        self.maxit = max_iterations
        self.tol = tolerance
        self.verbose = verbose

    def set_currents(self, j_cv, I):
        """
        Set the currents for the coils and vessel.
        
        @param j_cv: Current density in the vessel.
        @param I: List of currents for each coil.
        """
        if( len(I) != len(self.tags['coils'])): 
            raise ValueError("Number of currents does not match number of coils!")
        else:
            # Update the data in self.params:
            self.data['j_cv'] = j_cv
            self.data['j_coils'] = compute_j_coils(self.m, self.tags['coils'], I)


    def set_initial_guess(self, initial_guess=Constant(0.0)):
        """
        Set the initial guess for the flux function psi.
        
        @param initial_guess: Initial guess for psi, can be a constant or a function.
        """
        if isinstance(initial_guess, Function):
            self.psi.assign(initial_guess)
        elif isinstance(initial_guess, (int, float)):
            self.psi.interpolate(Constant(initial_guess))
        else:
            raise ValueError("Initial guess must be a constant or a Firedrake Function.")


    # CONSTRUCTOR METHODS:
    def import_mesh(self, geometry):
        """
        Set the mesh for the Grad-Shafranov problem, loaded from a .msh file in the 'meshes' folder.

        @param geometry: string identifying Tokamak geometry
        """

        path = "./meshes/" + geometry + ".msh"

        if not os.path.exists(path):
            raise FileNotFoundError(f"Mesh file {path} does not exist.")
        else: 
            print(f"\nLoading mesh for {geometry} geometry ...\n")
            self.m = Mesh(path, dim = 2, distribution_parameters={"partition": False}, reorder = True)
            self.m.init()


    # COMPUTE RESIDUAL:
    def compute_residual(self,psi_N):
        """
        Computes the residual of the Grad-Shafranov equation inside the plasma.
        This method calculates the residual based on the current flux function psi
        and the plasma parameters.

        @params psi_N: normalized flux as a Firedrake function.
        @return: The computed residual as a Firedrake Function.
        """
        # Apply Grad-Shafranov operator to psi:
        x,_ = SpatialCoordinate(self.m)
        mu0 = 4e-7 * pi 
        v = TestFunction(self.V)
        u = TrialFunction(self.V)
        #a = (1 / (mu0 * x) * dot(grad(u), grad(v))) * dx(domain=self.m)
        #A = assemble(a, bcs=self.dirichlet).petscmat

        Apsi = (1 / (mu0 * x) * dot(grad(self.psi), grad(v))) * dx(domain=self.m) - \
            self.plasma.domain_mask * self.data['G'](x, psi_N) * v * dx(self.tags['vacuum'], domain=self.m)
        res = Function(self.V)
        with assemble(Apsi, bcs=self.dirichlet).dat.vec_ro as v:
            res.dat.data[:] = v
      
        return sqrt(assemble(res**2 * dx(self.tags['vacuum'],domain=self.m)))


    # SOLVER METHODS:
    def normalize_flux(self, psi_N):
        """
        Normalize the poloidal flux function psi.
        The normalized flux ranges from 0 at the magnetic axis to 1 at the plasma boundary.

        @param psi_N: The normalized poloidal flux function to be updated.
        """
        # Normalize the poloidal flux:       
        denominator = self.plasma.psi_ma - self.plasma.psi0
        if denominator < 1e-14: # Avoid division by zero
            psi_N.assign(Constant(0.0))
        else:
            psi_N.interpolate((self.plasma.psi_ma - self.psi) / denominator)


    def solve(self):
        """
        Solve the Grad-Shafranov problem.
        """
        print("Initializing Grad-Shafranov problem...\n")

        # Define function to store flux and normalized flux from previous iteration:
        psi_old = Function(self.V)
        psi_N = Function(self.V)
        self.normalize_flux(psi_N)

        # Initialize solver:
        if self.algorithm == "Picard":
            solver = Picard(self.V,self.data,self.tags,self.dirichlet)
            print("Using fixed point iterations...\n")
        else:
            solver = Newton(self.V,self.data,self.tags,self.dirichlet)
            print("Using Newton's method...\n")

        # Neumann boundary datum computed with BEM
        q = Function(self.V)
        #q.assign(self.neumann.compute_datum(psi_old))

        # Set up the solver parameters
        it = 0
        converged = False
        err_vector = []
        err_vector.append(1)
        #res_vector = []

        # Initialize residual:
        #res_norm = self.compute_residual(psi_N)
        #res_vector.append(res_norm)

        # Start .pvd file to store the solution
        outfile = VTKFile("./results/flux.pvd")
        outfile.write(self.psi)     # write initial condition:

        while not converged and it < self.maxit:

            # Initialize flux from the previous iteration:
            psi_old.assign(self.psi)

            # Solve the variational problem
            if self.algorithm == "Picard":
                args = (psi_N, self.plasma.domain_mask,q)
            else:
                args = (
                    psi_old, psi_N, self.plasma.n, self.plasma.domain_mask, self.plasma.boundary_mask,
                    self.plasma.psi0, self.plasma.psi_ma, self.plasma.x0_idx, self.plasma.x1_idx, q
                )
            self.psi.assign(solver.perform_iteration(*args))
            
            # Update the plasma mask (smoothed):
            self.plasma.update(self.psi)

            # Normalize flux:
            self.normalize_flux(psi_N)

            # Compute residual:
            #res_norm = self.compute_residual(psi_N)
            #es_vector.append(res_norm)

            # Compute H1 error:
            self.err = errornorm(self.psi, psi_old, 'L2') / norm(psi_old, 'L2')
            err_vector.append(self.err)

            # Update Neumann BC datum:
            #q.assign(self.neumann.compute_datum(self.psi))

            # Print iteration information
            if self.verbose:
                print(f"Iteration {it+1}: H1 Error = {self.err:.6e}, psi at boundary = {self.plasma.psi0:.6f}, max psi = {self.plasma.psi_ma:.6f}")
            
            outfile.write(self.psi)
            # Check convergence:
            if self.err < self.tol:
                converged = True
            else:
                it += 1

        self.plot_convergence(err_vector)
        if not converged:
            print(f"Solver did not converge after {self.maxit} iterations. Final L2 norm of residual = {err_vector[-1]:.6e}\n")
        else:
            print(f"Solver converged in {it} iterations. Final L2 norm of residual = {err_vector[-1]:.6e}\n")


    # PLOT METHODS:
    def plot_mesh(self):
        """
        Display the mesh of the Grad-Shafranov problem.
        This method saves a plot of the mesh in the results directory.
        """
        self.plot.mesh()

    def plot_flux(self, path = "./results/flux_plot.png"):
        """
        Plot the flux function psi.
        This method saves a contour plot of the flux function in the results directory.

        @params path: path to save the plot
        """
        self.plot.flux(self.psi, self.plasma.psi0, self.plasma.d, self.plasma.x1_idx, self.plasma.x0_idx, path, self.plasma.psi_X_point)

    def plot_convergence(self, err_vector):
        """
        Plot the error behavior of the Grad-Shafranov equation.

        @params err_vector: vector containing the error/residual value at each iteration
        """
        path = "./results/error_plot.png"
        fig, ax = plt.subplots()
        ax.plot(range(1, len(err_vector)), err_vector[1:], marker='o')
        ax.set_yscale('log')
        ax.set_xlabel('Num. of iteration')
        ax.set_ylabel('Error')
        ax.set_title('Newton Iterations')
        ax.grid(True)
        plt.savefig(path)
        plt.close()

        path="./results/convergence_history.txt"
        with open(path, "w") as f:
            f.write("Error\n")
            np.savetxt(f, err_vector, fmt="%.6e")