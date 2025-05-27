from firedrake import *
from firedrake.pyplot import triplot, tricontourf
from .functions.geometry import generate_mesh
from .functions.regions import initialize_plasma_mask, define_vessel_mask, define_coils_mask, update_plasma_mask
from .functions.varf import GS_varf_Picard
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

        # Set boundary conditions
        self.BCs = [DirichletBC(self.V, 0.0, 'on_boundary')]

        # Identify the different regions in the mesh:
        self.plasma_mask = initialize_plasma_mask(self.params["vessel"], self.V, self.x, self.y)
        self.vessel_mask = define_vessel_mask(self.params["vessel"], self.params["j_cv"], self.V, self.x, self.y)
        self.coils_mask = define_coils_mask(self.params["I"], self.params["coils"], self.V, self.x, self.y)

        self.converged = False
        self.algorithm = "Picard" # Default algorithm for solving the Grad-Shafranov equation

    # SETTERS:
    def set_algorithm(self, algorithm):
        """
        Set the algorithm for solving the Grad-Shafranov equation.
        
        :param algorithm: The algorithm to use (e.g., "Picard", "Newton").
        """
        if algorithm not in ["Picard"]:
            raise ValueError("Invalid algorithm. Choose 'Picard' or 'Newton'.")
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
            # Update current density in vessel_mask:
            self.vessel_mask.assign(Constant(j_cv/self.params["j_cv"]) * self.vessel_mask)    

            # Update the currents in the coils_mask:
            self.coils_mask = define_coils_mask(I, self.params["coils"], self.V, self.x, self.y)

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
        geometry = self.params.get("geometry", "build")

        # Either build the mesh based on the geometry parameters provided
        if( geometry == "build"):
            path = "./meshes/generated_mesh.msh"
            
            geometry_params = {
                "x0": self.params.get("x0", 0.0),
                "y0": self.params.get("y0", 0.0),
                "R": self.params.get("R", 1.0),
                "thickness": self.params.get("thickness", 0.1),
                "coils": self.params.get("coils", []),
                "domain_size": self.params.get("domain_size", (1.0, 1.0)),
                "mesh_size_min": self.params.get("mesh_size_min", 0.01),
                "mesh_size_max": self.params.get("mesh_size_max", 0.05),
                "vacuum_mesh_size": self.params.get("vacuum_mesh_size", 0.02),
                "limiter_pts": self.params.get("limiter_pts", None),
            }

            # Generate the mesh using the geometry parameters
            generate_mesh(geometry_params, path)
            self.Mesh = Mesh(path)
            self.limiter = geometry_params["limiter_pts"]

        # Or load a pre-defined mesh from a file
        else:
            path = "./meshes" + geometry + ".msh"
            if not os.path.exists(path):
                raise FileNotFoundError(f"Mesh file {path} does not exist.")
            else: 
                self.Mesh = Mesh(path)

            # Set the limiter points for the plasma boundary:
            self.limiter = DirichletBC(V, 0.0, 'limiter').nodes          # case limiter is a closed line
    

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
        self.phi = TestFunction(self.V)  # Test function
        self.psi = Function(self.V)       # Trial flux function


    def assemble_system(self, phi, psi_old):
        """
        Assemble the system of equations for the Grad-Shafranov problem.
        This method should be implemented to set up the variational problem.
        """

        if self.algorithm == "Picard":
            a = GS_varf_Picard(self.x, self.y, self.params["G"], self.phi, self.psi, psi_old, 
                                    self.plasma_mask, self.vessel_mask, self.coils_mask)
            
        else:
            NotImplementedError("Method not implemented. Please choose a valid method.")

        return a

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

        # Assemble the system
        a = self.assemble_system(self.phi, psi_old)

        # Set up the solver parameters
        maxit = self.params.get("max_iterations", 100)
        tol = self.params.get("tolerance", 1e-6)
        it = 0
        self.converged = False

        # Start file to store the solution
        outfile = VTKFile("./results/flux.pvd")

        while not self.converged and it < maxit:
            # Solve the variational problem
            solve(a == 0, self.psi, bcs = self.BCs)

            # Compute error:
            self.err = errornorm(self.psi, psi_old, 'H1') / norm(psi_old, 'H1')
            
            # Update psi_old and mask for the next iteration
            psi_old.assign(self.psi)
            self.psi0 = update_plasma_mask(self.psi, self.limiter, self.plasma_mask)

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
        plt.scatter(*zip(*self.limiter), color='blue', label='Limiter Points')

        if( self.params.get("geometry", "build") == "build" ):
            # Plot coils 
            for coil in self.params["coils"]:
                x_min, x_max, y_min, y_max = coil
                plt.plot([x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min], color='red', label='Coil Edges')
            # Plot vessel
            Vessel = self.params["vessel"]
            inner_wall = plt.Circle((Vessel[0], Vessel[1]), Vessel[2], color='yellow', fill=False, label='Vessel')
            outer_wall = plt.Circle((Vessel[0], Vessel[1]), Vessel[2] + Vessel[3], color='yellow', fill=False, linestyle='dashed')
            ax.add_artist(inner_wall)
            ax.add_artist(outer_wall)

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
        #fig.colorbar(tripcolor(self.psi, axes=ax))
        tricontourf(self.psi, levels=20, axes=ax)
        plt.scatter(*zip(*self.limiter), color='blue', label='Limiter Points')
        # Include domain structures:
        if( self.params.get("geometry", "build") == "build" ):
            # Plot coils 
            for coil in self.params["coils"]:
                x_min, x_max, y_min, y_max = coil
                plt.plot([x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min], color='black', label='Coil Edges')
            # Plot vessel
            Vessel = self.params["vessel"]
            inner_wall = plt.Circle((Vessel[0], Vessel[1]), Vessel[2], color='red', fill=False, label='Vessel')
            outer_wall = plt.Circle((Vessel[0], Vessel[1]), Vessel[2] + Vessel[3], color='red', fill=False, linestyle='dashed')
            ax.add_artist(inner_wall)
            ax.add_artist(outer_wall)
        plt.title(r"Solution $\psi$. In red the plasma boundary")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.axis('equal')
        plt.savefig("./results/contour_plot.png")
        plt.close()