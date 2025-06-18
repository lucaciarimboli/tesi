from firedrake import *
from firedrake.pyplot import triplot, tricontour, tricontourf
from .functions.geometry import generate_mesh
from .functions.mesh_tags import get_tags
from .functions.varf import Picard_varf, Newton_varf
from .functions.coils import compute_j_coils
import os
import matplotlib.pyplot as plt
import numpy as np

class GradShafranovSolver:

    # CONSTRUCTOR:
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
        self.Mesh.init()
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
        
        # Compute the coils current density:
        self.j_coils = compute_j_coils(self.Mesh, self.tags['coils'], params["I"])

        self.converged = False
        self.set_algorithm(params.get("algorithm", "Picard"))  # Default algorithm is Picard


    # SETTERS:
    def set_algorithm(self, algorithm):
        """
        Set the algorithm for solving the Grad-Shafranov equation.

        :param algorithm: The algorithm to use (e.g., "Picard", "Newton").
        """
        if algorithm not in ["Picard", "Marder-Weitzner", "Newton"]:
            self.algorithm = "Picard"  # Default to Picard if invalid algorithm is provided
            print(f"Invalid algorithm '{algorithm}' provided. Defaulting to 'Picard'.")
        else:
            self.algorithm = algorithm

    def set_iterations_params(self, max_iterations = None, tolerance = None, verbose=False):
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
            self.j_coils = compute_j_coils(self.Mesh, self.tags['coils'], I)

    def set_initial_guess(self, initial_guess=Constant(0.0)):
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

            print("\nBuilding mesh based on custom geometry parameters...")
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
            print(f"\nLoading mesh for geometry '{self.params['geometry']}'...")
            path = "./meshes/" + self.params['geometry'] + ".msh"
            #if not os.path.exists(path):
            #    raise FileNotFoundError(f"Mesh file {path} does not exist.")
            #else: 
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
        print(f"Defining function spaces of '{family}{degree}' type...")
        self.V = FunctionSpace(self.Mesh, family, degree)
        self.x, self.y = SpatialCoordinate(self.Mesh)

        # Define test and trial functions:
        self.psi_trial = TrialFunction(self.V)  # Trial function
        self.phi = TestFunction(self.V)         # Test function

        self.psi = Function(self.V)             # Result flux function
        self.plasma_mask = Function(self.V)     # Plasma mask function


    # SOLVER METHODS:
    def normalize_flux(self, psi, psi_N):
        """
        Normalize the poloidal flux function psi.
        The normalized flux ranges from 0 at the magnetic axis to 1 at the plasma boundary.
        To avoid numerical instabilities the maximum value of the normalized flux is set to 0.99.
        Starting from the magnetic surface which correspond to a normalized psi of 0.99,
        the normalized flux is indeed set as constant = 0.99.

        :param psi: The poloidal flux function to normalize.
        :param psi_N: The normalized poloidal flux function to be updated.
        """
        # Compute the flux value in the magnetic axis:
        self.psi_max = psi.vector().max()

        # Compute the flux value at the plasma boundary:
        if isinstance(self.limiter[0], tuple):
            self.psi0 = max(psi.at(pt) for pt in self.limiter)
        else:
            self.psi0 = max(psi.dat.data[self.limiter])
 
        # Normalize the poloidal flux:       
        denominator = self.psi_max - self.psi0
        if denominator < 1e-14: # Avoid division by zero
            psi_N.assign(Constant(0.0))
            self.denom = 1
        else:
            psi_N.interpolate(conditional(
                (self.psi_max - psi) / denominator > 0.99,
                Constant(0.99),
                (self.psi_max - psi) / denominator
            )) 
            self.denom = denominator

    def perform_iteration(self, psi_N, psi_old, solver_params):
        '''
        Perform one iteration of the solution of Grad-Shafranov equation,
        using the iterative method prescribed by the class status.

        :param psi_N: normalized flux field at the previous iteration
        :param psi_old: magnetic flux field at the previous iteration
        :param solver_params: eventually contains other needed parameters
                (e.g. relaxation coeff for multi-step methods)
        '''

        # Case Picard:
        if self.algorithm == "Picard":
            a,L = Picard_varf(self.Mesh, self.x, self.params["G"], self.phi, self.psi_trial, psi_N,
                            self.plasma_mask, self.params["j_cv"], self.j_coils,
                            self.tags['vacuum'], self.tags['vessel'], self.tags['coils'])
            solve(a == L, self.psi, bcs = [self.BCs])
            
            # Per ottimizzare potrei assemblare "A" fuori dal ciclo while:
            '''
            a = Picard_Varf_a( etc. etc. )
            A = assemble(a, bcs=[self.BCs])
            solver = LinearSolver(A,
                  solver_parameters={'ksp_type': 'preonly',
                                     'pc_type': 'lu'})
            '''
            # E poi aggiornare ad ogni iter solo b:
            '''
            L = Picard_varf_L( etc. etc. )
            b = assemble(L, bcs=[self.BCs])
            solver.solve(self.psi, b)
            '''

        elif self.algorithm == "Marder-Weitzner":
            
            alpha = solver_params["alpha"]     # relaxation parameter:
            psi_step = Function(self.V) # function to store intermediate step solution:
            psi_step.assign(psi_old)

            # Compute first-step flux of Marder-Weitzner method:
            a,L = Picard_varf(self.Mesh, self.x, self.params["G"], self.phi, self.psi_trial, psi_N,
                            self.plasma_mask, self.params["j_cv"], self.j_coils,
                            self.tags['vacuum'], self.tags['vessel'], self.tags['coils'])
            solve(a == L, psi_step, bcs = [self.BCs])

            # Normalize the intermediate step flux
            self.normalize_flux(psi_step, psi_N)
            #self.plasma_mask.interpolate(0.5 + 0.5 * tanh((self.psi_step - self.psi0) / (epsilon * self.psi0)))

            # Compute second step flux of Marder-Weitzner method:
            a,L = Picard_varf(self.Mesh, self.x, self.params["G"], self.phi, self.psi_trial, psi_N,
                            self.plasma_mask, self.params["j_cv"], self.j_coils,
                            self.tags['vacuum'], self.tags['vessel'], self.tags['coils'])
            solve(a == L, self.psi, bcs = [self.BCs])

            # Update flux using the Marder-Weitzner method with relaxation alpha:
            self.psi.assign( (1-alpha) * psi_old + 2*alpha * psi_step - alpha * self.psi )

        elif self.algorithm == "Newton":
            print(f'\nValue of the denominator: {self.denom}')
            a, L = Newton_varf(self.Mesh, self.x, self.params["G"], self.phi, self.psi_trial, psi_N, psi_old,
                            self.denom, self.plasma_mask, self.params["j_cv"], self.j_coils,
                            self.tags['vacuum'], self.tags['vessel'], self.tags['coils'])
            #print("Max of residual vector F:", assemble(F).dat.data.max())
            #print("Norm of Jacobian matrix J:", assemble(J).petscmat.norm())
                
            solve(a == L, self.psi, bcs = [self.BCs])
            #self.psi.assign(self.psi + psi_old)

    def solve(self):
        """
        Solve the Grad-Shafranov problem.
        
        :return: The solution to the Grad-Shafranov equation.
        """
        print("Initializing Grad-Shafranov problem...\n")
        # Poloidal flux old iteration:
        psi_old = Function(self.V)      # poloidal flux old iteration
        psi_N = Function(self.V)        # normalized poloidal flux

        if "initial_guess" in self.params:
            psi_old.assign(self.params["initial_guess"])
            self.normalize_flux(psi_old, psi_N)
        else:
            print("No initial guess provided, using default values...")
            psi_old.interpolate(Constant(1e-5))  # Not 0 to avoid division by zero in the first iteration
            psi_N.interpolate(Constant(0.0))
            self.denom = 1

        # Initialize the plasma mask:
        #self.plasma_mask.interpolate(Constant(1.0))
        self.plasma_mask.interpolate(Constant(0.5))

        epsilon = 0.01  # Smoothing parameter for the plasma mask

        # Define function for intermediate solution of multi-step method
        solver_params = {}
        if self.algorithm == "Marder-Weitzner":
            solver_params["alpha"] = self.params.get("alpha", 0.5) # Relaxation parameter

        # Set up the solver parameters
        maxit = self.params.get("max_iterations", 100)
        tol = self.params.get("tolerance", 1e-6)
        it = 0
        self.converged = False

        # Start .pvd file to store the solution
        outfile = VTKFile("./results/flux.pvd")

        while not self.converged and it < maxit:

            # Solve the variational problem
            self.perform_iteration(psi_N, psi_old, solver_params)

            # Compute H1 error:
            self.err = errornorm(self.psi, psi_old, 'H1') / norm(psi_old, 'H1')
            
            # Normalize the poloidal flux:
            self.normalize_flux(self.psi, psi_N)

            # Update the plasma mask (smoothed:
            self.plasma_mask.interpolate(0.5 + 0.5 * tanh((self.psi - self.psi0) / (epsilon * self.psi0)))

            # Update psi_old and mask for the next iteration
            psi_old.assign(self.psi)

            # Print iteration information
            if(self.params.get("verbose", False)):
                print(f"Iteration {it+1}: H1 Error = {self.err:.6e}, psi at boundary = {self.psi0:.6f}, max psi = {self.psi_max:.6f}")

            # Write the solution to file
            outfile.write(self.psi)

            # Check convergence:
            if self.err < tol:
                self.converged = True
            else:
                it += 1

        if not self.converged:
            print(f"Solver did not converge after {maxit} iterations. Final H1 Error = {self.err:.6e}\n")
        else:
            print(f"Solver converged in {it} iterations. Final H1 Error = {self.err:.6e}\n")


    # PLOT METHODS:
    def display_mesh(self):
        """
        Display the mesh of the Grad-Shafranov problem.
        This method saves a plot of the mesh in the results directory.
        """
        print("Plotting mesh in file 'results/mesh_plot.png'...")
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
        This method saves a contour plot of the flux function in the results directory.
        """

        print("Plotting flux in file 'results/contour_plot.png'...")

        fig, ax = plt.subplots()
        tricontourf(self.psi, levels=50, cmap='viridis', axes=ax)
        tricontour(self.psi, levels=[self.psi0], colors='red', linewidths=2, axes=ax)
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


    # CONVERGENCE TEST:
    def compute_plasma_area(self):
        '''
        Computes the area of the plasma region in the poloidal plane.
        The computation uses Lebesgue measure, the integral is performed over the elements that are completely included by the plasma.
        '''
        psi0_func = Function(self.V)
        psi0_func.assign(Constant(self.psi0))
        plasma_region = SubDomainData(self.psi > psi0_func)

        return assemble(Constant(1.0) * dx(plasma_region,domain=mesh))

    def convergence_test(self, refinment_levels = 2):
        '''
        Performs a convercence test for the current set up.

        Other meshes are defined by unformly refining the solver's one at the current status.
        A simulation is performed for every mesh and the behavior of relevant quantities against the refinment level is plotted.

        The plot is saved in "results/grid_independence.png".
        The relevant quantites under inspection are:
            - (normalized) magnetic flux value at the magnetic axes;
            - (normalized) magnetic flux value at the plasma boundary;
            - (normalized) plasma region area in the poloidal plane.

        :param refinement_levels: number of mesh generated OTHER THAN the status one (default value is 2).
        '''

        psi_max_vect = []
        psi0_vect = []
        #plasma_area = []
        h_min = []

        print('Building meshes with increasing refinement levels...')
        meshes = MeshHierarchy(self.Mesh, refinment_levels)

        for n in range(refinment_levels):
            self.Mesh = meshes[n]
            self.Mesh.init()

            h_min.append(self.Mesh.cell_sizes.dat.data.min())
            print(f'\nSolving with minimum cell size h = {h_min[n]}...')

            # Adapt function spaces to the new refinment level:
            self.function_spaces()
            self.BCs = DirichletBC(self.V, 0.0, self.tags['boundary'])
            self.limiter = DirichletBC(self.V, 0.0, self.tags['limiter']).nodes

            self.solve()

            print('Saving results from simulation...')
            psi_max_vect.append(self.psi_max)
            psi0_vect.append(self.psi0)
            #plasma_area.append(self.compute_plasma_area())

        # Non-dimensionalize reference quantities:
        print(f'Plotting results...')
        h_min = np.array(h_min)
        #plasma_area = np.array(plasma_area) / plasma_area[0]
        psi_max_vect = np.array(psi_max_vect) / psi_max_vect[0]
        psi0_vect = np.array(psi0_vect) / psi0_vect[0]

        plt.figure()
        plt.plot(h_min, psi_max_vect, marker='o', label='psi_max')
        plt.plot(h_min, psi0_vect, marker='s', label='psi0')
        #plt.plot(h_min, plasma_area, marker='^', label='Relative Plasma Area')
        plt.xlabel('Minimum cell size h')
        plt.ylabel('Value')
        plt.title('Grid Independence Study')
        plt.legend()
        plt.grid(True)
        plt.savefig("./results/grid_independence.png")