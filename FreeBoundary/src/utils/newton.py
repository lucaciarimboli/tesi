from firedrake import *
from src.utils.functions.mask import delta_point

class Newton:
    '''
        Class to perform one fixed point iteration of the Grad-Shafranov problem.
    '''

    def __init__(self, V, data, tags, bcs):

        # Save data and mesh tags:
        self.j_plasma = data['G']
        self.j_plasma_derivative = data['dGdpsiN']
        self.j_coils = data['j_coils']
        self.j_cv = data['j_cv']
        self.tags = tags

        # Define mixed function space
        m = V.mesh()
        x, y = SpatialCoordinate(m)
        R = FunctionSpace(m, 'R', 0)
        self.W = V * R * R

        self.psi, self.lambda0, self.lambda1 = TrialFunctions(self.W)
        self.phi, self.v0, self.v1 = TestFunctions(self.W)

        mu0 = 4e-7 * pi  # Permeability of free space (in IS units)
        self.a = ( 1 / (mu0 * x) * dot(grad(self.psi), grad(self.phi)) ) * dx(domain=m)

        # Boundary conditions:
        self.bcs = DirichletBC(self.W.sub(0), bcs.function_arg, bcs.sub_domain)
        #self.bcs = DirichletBC(self.W.sub(0), 0.0, tags['boundary'])


    def linear_form(self, psi_N, plasma_mask):

        m = self.W.mesh()
        x,y = SpatialCoordinate(m)

        # Plasma current density contribution:
        L = plasma_mask *  self.j_plasma(x,psi_N) * self.phi * dx(self.tags['vacuum'], domain=m)

        # Vessel wall current density contribution:
        L += self.j_cv * self.phi * dx(self.tags['vessel'], domain=m)

        # Coils current density contribution:
        coils_tags = self.tags['coils']
        for i in range(len(self.j_coils)):
            L += self.j_coils[i] * self.phi * dx(coils_tags[i], domain=m)

        return L
    

    def Jacobian_form(self, psi_old, psi_N, dpsidn, plasma_mask, boundary_mask, psi0, psi_ax, x0_idx, x1_idx):

        m = self.W.mesh()
        x,y = SpatialCoordinate(m)

        # Gateaux derivative of psi_N w.r.t. psi:
        e = 1/(psi_ax - psi0)**2 * ((psi0-psi_ax)*self.psi + (psi_old-psi0)*self.lambda1 + (psi_ax-psi_old)*self.lambda0)

        # Surface integral term:
        d = plasma_mask * self.j_plasma_derivative(x,psi_N) * e * self.phi * dx(self.tags['vacuum'], domain = m)

        # Boundary term:
        d += boundary_mask * self.j_plasma(x,1.0) * 1/dpsidn * (self.lambda0 - self.psi) * self.phi * dx(self.tags['vacuum'], domain = m)

        # Point sources in x0 and x1 for trial function evaluation:
        V = self.W.sub(0)
        delta_x0 = delta_point(V,x0_idx)
        delta_x1 = delta_point(V,x1_idx)

        # Bilinear forms for the lambda values:
        l0 = ( self.psi - self.lambda0 ) * delta_x0 * self.v0 * dx(domain = m, degree=20)
        l1 = ( self.psi - self.lambda1 ) * delta_x1 * self.v1 * dx(domain = m, degree=20)
        # (degree = 20 is need to avoid the warning about interpolation)

        return - d + l0 + l1


    def perform_iteration(self, psi_old, psi_N, n, plasma_mask, boundary_mask, psi0, psi_ax, x0_idx, x1_idx):

        # Normal derivative for psi at plasma boundary (computed everywhere)
        V = self.W.sub(0)
        dpsidn = Function(V).interpolate(inner(grad(psi_old),n))

        # Define linear and bilinear forms:
        a = self.a + self.Jacobian_form(psi_old, psi_N, dpsidn, plasma_mask, boundary_mask, psi0, psi_ax, x0_idx, x1_idx)
        L = self.linear_form(psi_N, plasma_mask)

        # Solve for psi:
        w = Function(self.W)
        solve(a == L, w, bcs = self.bcs)

        return w.sub(0)