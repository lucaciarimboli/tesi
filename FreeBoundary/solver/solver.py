from firedrake import *
from firedrake.pyplot import triplot, tricontour, tricontourf
from .functions.geometry import generate_mesh
from .functions.mesh_tags import get_tags
from .functions.varf import Picard_varf
import os
import matplotlib.pyplot as plt

class GradShafranovSolver:

    def __init__(self, params):
        """
        Constructor for the GradShafranovSolver class.
        This class is designed to solve the Grad-Shafranov equation using a finite element method.
        It initializes the solver with the given parameters and prepares the mesh for the problem.
        
        :param params: Dictionary containing the parameters for the solver.
        """
        self.params = params

        # Initialize the mesh and function spaces
        self.build_mesh()
        self.function_spaces()

        # Identify the tags for the mesh:
        geometry = params.get("geometry", "custom")
        if geometry == "custom":
            self.tags = get_tags(geometry, params)
        else:
            self.tags = get_tags(geometry)

        # Set boundary conditions
        self.BCs = DirichletBC(self.V, 0.0, self.tags['boundary'])

        # Extract the nodes that lie on the limiter:

        if self.tags['limiter'] is not None:
            self.limiter = DirichletBC(self.V, 0.0, self.tags['limiter']).nodes

        if self.tags['limiter_pts'] is not None:
            if params.get("limiter_pts", None) is not None:
                self.limiter = params["limiter_pts"]
        # ------------------------------------------------------------------------------------------------------------------------- #
        # AGGIUNGERE UN METODO PER "ESTRARRE" I SINGOLI PHYSICAL POINTS USANDO I TAG, SENZA DUNQUE PASSARE LE COORDINATE DIRETTAMENTE
        # ------------------------------------------------------------------------------------------------------------------------- #
        
        self.converged = False
        self.set_algorithm(params.get("algorithm", "Picard"))  # Default algorithm is Picard


    # SETTERS:
    def set_algorithm(self, algorithm):
        """
        Set the algorithm for solving the Grad-Shafranov equation.

        :param algorithm: The algorithm to use (e.g., "Picard", "Marder-Weitzner").
        """
        if algorithm not in ["Picard", "Marder-Weitzner"]:
            self.algorithm = "Picard"  # Default to Picard if invalid algorithm is provided
            print(f"Invalid algorithm '{algorithm}' provided. Defaulting to 'Picard'.")
        else:
            self.algorithm = algorithm

    def set_iterations_params(self, max_iterations, tolerance, verbose=False):
        """
        Set the parameters for the solver iterations.
        
        :param max_iterations: Maximum number of iterations for the solver.
        :param tolerance: Tolerance for convergence.
        :param verbose: If True, print iteration details.
        """
        self.params["max_iterations"] = max_iterations
        self.params["tolerance"] = tolerance
        self.params["verbose"] = verbose

    def set_currents(self, j_cv, I):
        """
        Set the currents for the coils and vessel.
        
        :param j_cv: Current density in the vessel.
        :param I: List of currents for each coil.
        """
        if( len(I) != len(self.params["coils"]) ): 
            raise ValueError("Number of currents does not match number of coils!")
        else:
            # Update the data in self.params:
            self.params["j_cv"] = j_cv
            self.params["I"] = I


    def set_initial_guess(self, initial_guess):
        """
        Set the initial guess for the flux function psi.
        
        :param initial_guess: Initial guess for psi, can be a constant or a function.
        """
        if isinstance(initial_guess, Function):
            self.params["initial_guess"] = initial_guess
        elif isinstance(initial_guess, (int, float)):
            self.params["initial_guess"] = Constant(initial_guess)
        else:
            raise ValueError("Initial guess must be a constant or a Firedrake Function.")


    # METHODS FOR THE CONSTRUCTOR:
    def build_mesh(self):
        """
        Set the mesh for the Grad-Shafranov problem.
        The mesh can be either loaded from the provided ones in the folder "meshes",
        or can be generated specifying the geometry parameters. 
        """

        # Either build the mesh based on the geometry parameters provided
        if( self.params.get("geometry", "custom") == "custom" ):

            path = "./meshes/custom_tokamak.msh"
            
            geometry_params = {
                "boundary": self.params.get("boundary", [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]), 
                "outer_wall": self.params.get("vessel_outer_wall", None),
                "inner_wall": self.params.get("vessel_inner_wall", None),
                "coils": self.params.get("coils", []),
                "limiter_pts": self.params.get("limiter_pts", None),
                "limiter_line": self.params.get("limiter_line", None),
                "mesh_size_min": self.params.get("mesh_size_min", 0.01),
                "mesh_size_max": self.params.get("mesh_size_max", 0.05),
                "limiter_mesh_size": self.params.get("limiter_mesh_size", 0.02),
                "limiter_dist_max": self.params.get("dist_from_limiter", 0.0),
                "coils_mesh_size": self.params.get("coils_mesh_size", self.params.get("mesh_size_min")),
                "coils_dist_max": self.params.get("dist_from_coils", 0.0),
            }

            # Generate the mesh using the geometry parameters
            generate_mesh(geometry_params, path)
            self.Mesh = Mesh(path, distribution_parameters={"partition": False})

        # Or load a pre-defined mesh from a file
        else:
            path = "./meshes/" + self.params['geometry'] + ".msh"
            if not os.path.exists(path):
                raise FileNotFoundError(f"Mesh file {path} does not exist.")
            else: 
                self.Mesh = Mesh(path, distribution_parameters={"partition": False})
    

    def function_spaces(self, family=None, degree=None):
        """
        Define the function spaces for the Grad-Shafranov problem.

        Future edit: if adding dependcy on mu, need to define more spaces!        
        """
        # Define the function space (CG by default)
        if family is None:
            family = self.params.get("function_space_family", "CG")
        if degree is None:
            degree = self.params.get("function_space_degree", 1)
        self.V = FunctionSpace(self.Mesh, family, degree)
        self.x, self.y = SpatialCoordinate(self.Mesh)

        # Define test and trial functions:
        self.phi = TestFunction(self.V)         # Test function
        self.psi = Function(self.V)             # Trial flux function
        self.plasma_mask = Function(self.V)     # Plasma mask function


    # SOLVER METHOD:
    def solve(self):
        """
        Solve the Grad-Shafranov problem.
        
        :return: The solution to the Grad-Shafranov equation.
        """
        
        # Set initial guess for psi
        psi_old = Function(self.V)  
        if "initial_guess" in self.params:
            psi_old.interpolate(self.params["initial_guess"])
        else:
            psi_old.interpolate(Constant(0.0))

        # Initialize the plasma mask:
        self.plasma_mask.interpolate(Constant(1.0))
        epsilon = 0.01  # Smoothing parameter for the plasma mask

        # Define intermidiate solutions for multi-step methods if needed:
        if self.algorithm == "Marder-Weitzner":
            psi_step = Function(self.V)
            alpha = self.params.get("alpha", 0.5)  # Default alpha value for Marder-Weitzner method

        # Set up the solver parameters
        maxit = self.params.get("max_iterations", 100)
        tol = self.params.get("tolerance", 1e-6)
        it = 0
        self.converged = False

        # Start file to store the solution
        outfile = VTKFile("./results/flux.pvd")

        while not self.converged and it < maxit:

            # Solve the variational problem

            if self.algorithm == "Picard":
                a = Picard_varf(self.Mesh, self.x, self.params["G"], self.phi, self.psi, psi_old,
                                self.plasma_mask, self.params["j_cv"], self.params["I"],
                                self.tags['vacuum'], self.tags['vessel'], self.tags['coils'])
                solve(a == 0, self.psi, bcs = [self.BCs])

            elif self.algorithm == "Marder-Weitzner":
                a = Picard_varf(self.Mesh, self.x, self.params["G"], self.phi, psi_step, psi_old,
                                 self.plasma_mask, self.params["j_cv"], self.params["I"],
                                 self.tags['vacuum'], self.tags['vessel'], self.tags['coils'])
                solve(a == 0, psi_step, bcs = [self.BCs])
                a = Picard_varf(self.Mesh, self.x, self.params["G"], self.phi, self.psi, psi_step,
                                 self.plasma_mask, self.params["j_cv"], self.params["I"],
                                 self.tags['vacuum'], self.tags['vessel'], self.tags['coils'])
                solve(a == 0, self.psi, bcs = [self.BCs])
                self.psi.assign( (1-alpha) * psi_old + 2*alpha *psi_step - alpha * self.psi )

            else:
                raise ValueError(f"Unknown algorithm '{self.algorithm}'. Supported algorithms are 'Picard' and 'Marder-Weitzner'.")

            # Compute error:
            self.err = errornorm(self.psi, psi_old, 'H1') / norm(psi_old, 'H1')
            
            # Update psi_old and mask for the next iteration
            psi_old.assign(self.psi)
            
            # Update the plasma mask based on the new psi
            if isinstance(self.limiter[0], tuple):  # List of (x, y) coordinates
                self.psi0 = max(self.psi.at(pt) for pt in self.limiter)
            else:  # List/array of node indices
                self.psi0 = max(self.psi.dat.data[self.limiter])
            self.plasma_mask.interpolate(0.5 + 0.5 * tanh((self.psi - self.psi0) / (epsilon * self.psi0)))

            # Print iteration information
            if(self.params.get("verbose", False)):
                print(f"Iteration {it+1}: H1 Error = {self.err:.6e}, psi at boundary = {self.psi0:.6f}, max psi = {self.psi.vector().max():.6f}")

            # Write the solution to file
            outfile.write(self.psi)

            # Check convergence:
            if self.err < tol:
                self.converged = True
            else:
                it += 1

        if not self.converged:
            print(f"Solver did not converge after {maxit} iterations. Final H1 Error = {self.err:.6e}")
        else:
            print(f"Solver converged in {it} iterations. Final H1 Error = {self.err:.6e}")

    # PLOT METHODS:
    def display_mesh(self):
        """
        Display the mesh of the Grad-Shafranov problem.
        This method can be implemented to visualize the mesh.
        """

        fig, ax = plt.subplots()
        triplot(self.Mesh, axes=ax)
        #plt.scatter(*zip(*self.limiter), color='blue', label='Limiter Points')

        #if( self.params.get("geometry", "build") == "build" ):
            # Plot coils 
        #    for coil in self.params["coils_adapted_format"]:
        #        x_min, x_max, y_min, y_max = coil
        #        plt.plot([x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min], color='red', label='Coil Edges')
            # Plot vessel
        #    Vessel = self.params["vessel"]
        #    inner_wall = plt.Circle((Vessel[0], Vessel[1]), Vessel[2], color='yellow', fill=False, label='Vessel')
        #    outer_wall = plt.Circle((Vessel[0], Vessel[1]), Vessel[2] + Vessel[3], color='yellow', fill=False, linestyle='dashed')
        #    ax.add_artist(inner_wall)
        #    ax.add_artist(outer_wall)

        plt.title(r"Domain")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis('equal')
        plt.savefig("./results/mesh_plot.png")
        plt.close()

    def plot_flux(self):
        """
        Plot the flux function psi.
        This method can be implemented to visualize the solution.
        """

        fig, ax = plt.subplots()
        tricontourf(self.psi, levels=50, cmap='viridis', axes=ax, figure=fig)
        tricontour(self.psi, levels=[self.psi0], colors='red', linewidths=2, figure=fig, axes=ax)
        #triplot(self.Mesh, axes=ax)
        # Plot the plasma boundary:
        #plt.contour(self.x, self.y, self.psi, levels=[self.psi0], colors='red')
        if self.params.get("limiter_pts", None) is not None:
            # self.limiter is a list of (x, y) coordinates
            if isinstance(self.limiter, (list, tuple)) and isinstance(self.limiter[0], tuple):
                plt.scatter(*zip(*self.limiter), color='blue', label='Limiter Points')
            else:
                # self.limiter is a list/array of node indices: plot their coordinates
                coords = self.V.mesh().coordinates.dat.data[self.limiter]
                plt.scatter(coords[:, 0], coords[:, 1], color='blue', label='Limiter Points')
        elif self.params.get("limiter", None) is not None:
            # self.limiter is a list/array of node indices: plot their coordinates
            coords = self.V.mesh().coordinates.dat.data[self.limiter]
            plt.plot(coords[:, 0], coords[:, 1], '-', color='blue', label='Limiter')

        # Include domain structures:
        # [...]

        plt.title(r"Solution $\psi$. In red the plasma boundary")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.axis('equal')
        plt.savefig("./results/contour_plot.png")
        plt.close()