from firedrake import *

class Picard:
    '''
        Class to perform one fixed point iteration of the Grad-Shafranov problem.
    '''

    def __init__(self, V, data, tags, bcs):

        # Save data and mesh tags:
        self.j_plasma = data['G']
        self.j_coils = data['j_coils']
        self.j_cv = data['j_cv']
        self.tags = tags

        # Dirichlet boundary conditions:
        self.bcs = bcs

        # Define trial and test function:
        psi = TrialFunction(V)
        self.phi = TestFunction(V)

        # Define the bilinear form:
        m = V.mesh()
        x,_ = SpatialCoordinate(m)
        mu0 = 4e-7 * pi  # Permeability of free space (in IS units)
        a = ( 1 / (mu0 * x) * dot(grad(psi), grad(self.phi)) ) * dx(domain=m)

        # Define the linear solver:
        A = assemble(a,bcs=bcs)
        self.linear_solver = LinearSolver(A,solver_parameters={'ksp_type': 'preonly', 'pc_type': 'lu'})

    
    def linear_form(self, psi_N, plasma_mask, g):

        V = psi_N.function_space()
        m = V.mesh()
        x,_ = SpatialCoordinate(m)

        # Plasma current density contribution:
        L = plasma_mask *  self.j_plasma(x,psi_N) * self.phi * dx(self.tags['vacuum'], domain=m)

        # Vessel wall current density contribution:
        L += self.j_cv * self.phi * dx(self.tags['vessel'], domain=m)

        # Neumann bondary condition contribution:
        L += g * self.phi * ds(self.tags['neumann_boundary'], domain=m)

        # Coils current density contribution:
        coils_tags = self.tags['coils']
        for i in range(len(self.j_coils)):
            L += self.j_coils[i] * self.phi * dx(coils_tags[i], domain=m)

        return L
    

    def perform_iteration(self, psi_N, plasma_mask, g):
        
        # Update right hand side
        L = self.linear_form(psi_N, plasma_mask, g)
        b = assemble(L,bcs = self.bcs)

        # Solve for psi:
        V = psi_N.function_space()
        psi = Function(V)
        self.linear_solver.solve(psi,b)

        return psi