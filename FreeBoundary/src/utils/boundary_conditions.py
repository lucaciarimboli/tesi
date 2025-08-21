from firedrake import *
from firedrake.mesh import plex_from_cell_list
import numpy as np
from src.utils.functions.boundary_integrals import *


class JN_coupling_BCs:
    '''
        Solver for the boundary integral equation of Johnson-Nédélec.
        The solution of this equation is used as boundary condition on the artificial domain
        for the free-boundary Grad-Shafranov equation. 
    '''

    def __init__(self, V, boundary_tag):
        '''
            Constructor for the BEM solver class:
            @param boundary_tag: mesh tag for the artificial boundary
            @param psi: initial condition for the flux function
        '''

        self.boundary_tag = boundary_tag

        # Create 1D mesh on the boundary:
        m2d = V.mesh()
        self.boundary_mesh(m2d)
        self.m.init()
        
        # Function space on the boundary mesh:
        V_info = V.ufl_element()
        self.Q = FunctionSpace(self.m, V_info.family(), V_info.degree())

        # Initialize the integral function K:
        self.K_func = Function(self.Q)

        # Define boundary mask F and map "neighbor_segments",
        # which are needed to fill "K_func".
        self.integral_mask(m2d)

        # Compute Green functions and matrix M:
        self.assemble_Green_function(V)
        self.assemble_matrix(V)


## CONSTRUCTOR METHODS:
    def boundary_mesh(self, m2d):
        '''
            "Extract" a one-dimensional mesh for the boundary of the given two-dimensional mesh.
            Nodes indexing is preserved, external facets indexing is not.
            @param m2d: Two-dimensional mesh
        '''

        # Extract dofs coordinates:
        coord_func = Function(VectorFunctionSpace(m2d, "CG", 1)).interpolate(as_vector(SpatialCoordinate(m2d)))
        self.dof_coords = coord_func.dat.data_ro[:]

        # Fill "Q_dof_coords" vector with the coordinates of the boundary nodes of the 2D mesh
        self.Q_dofs = DirichletBC(FunctionSpace(m2d, "CG", 1), 0.0, self.boundary_tag).nodes
        Q_dof_coords = self.dof_coords[self.Q_dofs]

        # Identify extremes of the neumann boundary:
        self.extremes = np.argsort(Q_dof_coords[:, 0])[:2]

        # Fill "segments" with the 1D cells (segments) on the boundary
        n = len(self.Q_dofs)
        segments = set()   # format "set" to ignore 
        self.neighbors_map = []
        self.neighbors_dist = []

        for i in range(n):
            # Compute the distance with the others dofs:
            dist = np.zeros(n)
            for j in range(n):
                dist[j] = np.linalg.norm(Q_dof_coords[i]-Q_dof_coords[j])
            dist[i] = np.inf    # set distance with itself = inf

            # Identify the indexes of the two neighbouring dofs of dof_coords[i]
            # (one if the dof is an extreme)
            if i in self.extremes:
                neighbor_dofs = np.argmin(dist)
                segments.add(tuple(sorted([i,neighbor_dofs])))
                # Add neighbor idx and dist to maps as iterables (np.array)
                self.neighbors_map.append(np.array([self.Q_dofs[neighbor_dofs]]))
                self.neighbors_dist.append(np.array([dist[neighbor_dofs]]))

            else:    
                neighbor_dofs = np.argsort(dist)[:2]
                for k in neighbor_dofs:
                    segments.add(tuple(sorted([i,k])))
                # Add neighbors idexes and distances to maps
                self.neighbors_map.append(self.Q_dofs[neighbor_dofs])
                self.neighbors_dist.append(dist[neighbor_dofs])

        # Convert "segments" to array:
        segments = np.array(list(segments))

        # Build 1D boundary mesh from plex:
        plex = plex_from_cell_list(1, segments, Q_dof_coords, comm=m2d.comm)
        self.m = Mesh(plex, dim=1, reorder = False)


    def dofs_indexes(self, V):
        '''
            This function is needed in the constructor.
            It extracts and stores in "V_dofs" the indexes of the dofs that lie on the boundary triangles and on the coils.
            Over these nodes the green function will be evaluated in the computation of the integral functions K and L.
            This operation is done in order to avoid computing the green function over the whole domain since it is only needed
            at the boundary and in the coils + its gradient is needed at the boundary.

            @TODO: considera solo le boundary cells del neumann boundary
        '''
        # Extract the dof of the boundary mesh cells:
        m2d = V.mesh()
        bdry_cells = m2d.exterior_facets.facet_cell_map.values    #idx of boundary cells
        V_dofs = set()
        cell_nodes = V.cell_node_map()
        for dofs in cell_nodes.values[bdry_cells]:      # nodes format of kind: [[i,j,k]]
            for dof in dofs[0]:
                V_dofs.add(dof)

        # Per qualche motivo, non tutti i dof di "Q_dofs" sono inclusi in "V_dofs".
        # Per metterci una pezza -> faccio "a mano" l'aggiunta dei dof mancanti:
        for dof in self.Q_dofs:
            V_dofs.add(dof)

        return np.array(list(V_dofs))


    def assemble_Green_function(self,V):
        '''
            Assemble the Green functions on the boundary elements only.
            For each boundary dof, the list "G_list" contains at the corresponding
            index (considering 1d mesh indexing) with the fundamental solution
            obtained placing a point source in the dof itself
        '''

        # Identify dofs of the boundary cells:
        V_dofs = self.dofs_indexes(V)

        # List containing the Green function for each position of the point source 
        self.G_list = []

        for i in self.Q_dofs:
            X = self.dof_coords[i]

            # Fixed x, define G(x,y) on the boundary elements only
            G = Function(V) 
            for j in V_dofs:
                Y = self.dof_coords[j]
                G.dat.data[j] = Green_function(X[0],X[1],Y[0],Y[1])
            self.G_list.append(G)

            #-----------------------------------------------------------#
            #----- FINCHE' NON RISOLVO PROBLEMA DELL'INTEGRAZIONE-------#
            #-----------------------------------------------------------#
            G.dat.data[i] = 0.0
            #-----------------------------------------------------------#
            #-----------------------------------------------------------#


    def assemble_matrix(self,V):
        '''
            Assemble the Green function and the matrix M that allows to compute g from V(g).
        '''

        # Extract boundary dofs indexes (Q dofs indexes)
        mu0 = 4e-7 * pi

        # Matrix M for the computation of q:
        n = len(self.Q_dofs)
        self.M = np.zeros((n,n))

        for i in range(n):
            G = self.G_list[i]
            source_neighbors = self.neighbors_map[i]
            dof_i = self.Q_dofs[i]
            ri = self.dof_coords[dof_i][0]

            for j in range(n):
                dof_j = self.Q_dofs[j]
                rj = self.dof_coords[dof_j][0]
                Gj = G.dat.data_ro[dof_j]
                h = self.neighbors_dist[j]

                #self.M[i,j] = 0.0 # placeholder   
                # Point source dof on the diagonal (with singularity)
                if j==i:
                    if j in self.extremes:
                        self.M[i,j] = 1/mu0 * matrix_diagonal(h[0],ri)
                    else:
                        self.M[i,j] = 1/mu0 * ( matrix_diagonal(h[0],ri) + matrix_diagonal(h[1],ri) )

                # Handle singularity in neighboring dofs
                elif dof_j in source_neighbors:
                
                    if j in self.extremes:
                        self.M[i,j] = 1/mu0 * matrix_close(h[0],ri)
                    else:
                        # Re-order h so that the first element correspond to the singularity
                        j_neighbors = self.neighbors_map[j]
                        for k in range(2):
                            if j_neighbors[k] == dof_i:
                                h_close = h[k]
                            else:
                                h_far = h[k]

                        self.M[i,j] = 1/mu0 * (matrix_close(h_close,ri) + matrix_far(h_far,Gj,rj))

                # Simple trapezoidal quadrature far from singularity
                elif j in self.extremes:
                    self.M[i,j] = 1/mu0 * matrix_far(h[0],Gj,rj)
                else:
                    self.M[i,j] = 1/mu0 * (matrix_far(h[0],Gj,rj) + matrix_far(h[1],Gj,rj))      


    def integral_mask(self,m2d):
        '''
            Define a vector function F which is constant on boundary elements.
            The function is zero everywhere but in the mesh boundary, where is s.t.
            its outward normal is always =1.

            @param m: mesh 
 
            Inside boundary integrals with "assemble" one should include "*dot(F,n)", where n is the
            FacetNormal. By manually setting F = [0,0] on a boundary element it is possible
            to exclude selected boundary elements from the integration domain.

            To allow this operation, a map "nighbor_segments" is built. It provides for each boundary node
            of the mesh, the F dofs that lie on the two neighboring segments.

            @NOTE: The creation of an actual custom measure on the boundary that excludes given
            boundary elements based on their mesh indexes or on geometric considerations is not
            currently doable by means of Firedrake functions or by manipulating the plex.
            The introduction of this F might not be "elegant", but it is an effective and practical
            solution for excluding dofs with singularities from boundary integrals.
        '''
        W = VectorFunctionSpace(m2d, "HDivT", 0)
        self.F = Function(W)

        # Identify dofs that lie on boundary segments:
        boundary_seg_dofs = DirichletBC(FunctionSpace(m2d, "HDivT", 0),0.0,"on_boundary").nodes
        coord_func = Function(W).interpolate(as_vector(SpatialCoordinate(m2d)))
        boundary_coords = coord_func.dat.data_ro[boundary_seg_dofs]

        min_r = min(boundary_coords[:,0])
        max_r = max(boundary_coords[:, 0])
        max_z = max(boundary_coords[:, 1])
        min_z = min(boundary_coords[:, 1])

        for idx in boundary_seg_dofs:
            X = coord_func.dat.data_ro[idx]
            if X[1] == min_z:
                self.F.dat.data[idx][1] = -1
            if X[1] == max_z:
                self.F.dat.data[idx][1] = 1
            if X[0] == min_r:
                self.F.dat.data[idx][0] = -1
            if X[0] == max_r:
                self.F.dat.data[idx][0] = 1

        # Build map that associate to each Q dof, the two neighboring W dofs in the boundary
        n_dofs = len(self.Q_dofs)
        n_segments = len(boundary_seg_dofs)
        self.neighbor_segments = []

        for i in range(n_dofs):
            dof = self.Q_dofs[i]
            coords = self.dof_coords[dof]
            dist = np.zeros(n_segments)
            for j in range(n_segments):
                dist[j] = np.linalg.norm(coords-boundary_coords[j])
            ns_idxs = np.argsort(dist)[:2]
            self.neighbor_segments.append(boundary_seg_dofs[ns_idxs])


## METHODS FOR COMPUTING NEUMANN DATUM:
    def compute_K(self, psi):
        '''
            Update integral function K based on the value of the magnetic flux psi.
        '''

        V = psi.function_space()
        m2d = V.mesh()
        x,_ = SpatialCoordinate(m2d)
        n = FacetNormal(m2d)

        mu0 = 4e-7 * pi
        n_dofs = len(self.Q_dofs)

        # Assemble the new integral function K:
        for i in range(n_dofs):
            G = self.G_list[i]
            dof_i = self.Q_dofs[i]
            ns = self.neighbor_segments[i]

            # Mask-out neighbor segments of the i-esim boundary dof
            original_F_values = self.F.dat.data_ro[ns]
            self.F.dat.data[ns] = [[0.0,0.0],[0.0,0.0]]

            # Compute K(dof_i) contribution from all segments far from the singularity:
            integral = assemble( - 1/(mu0*x) * dot(grad(G),n) * dot(self.F,n) * psi * ds(self.boundary_tag,domain=m2d) )
            self.F.dat.data[ns] = original_F_values
            
            neighbors = self.neighbors_map[i]
            h = self.neighbors_dist[i]
            source_coords = self.dof_coords[dof_i]
            psi_source = psi.dat.data_ro[dof_i]
            
            # Add the contribution of the neighborhood of the singularity:
            for j, neigh in enumerate(neighbors):
                neigh_coords = self.dof_coords[neigh]
                if source_coords[0] == neigh_coords[0]:
                    psi_neigh = psi.dat.data_ro[neigh]
                    integral += 1/mu0 * K_neighborhood_integral(h[j],source_coords[0],psi_source,psi_neigh)

            # Update value of integral function K:
            self.K_func.dat.data[i] = integral 


    def solve_boundary_integral_eq(self,psi):
        '''
            Solves the boundary integral equation for V(q), solves the linear system to compute q from V(q).
        '''

        # Extract the trace of psi:
        psi_trace = Function(self.Q)
        psi_trace.dat.data[:] = psi.dat.data_ro[self.Q_dofs]
        
        # Define v(g) integral function
        v = Function(self.Q).assign(0.5 * psi_trace + self.K_func)

        # Solve the linear system to obtain g from V(q):
        return np.linalg.solve(self.M, v.dat.data_ro[:])
    

    def compute_datum(self,psi):
        '''
            Given the flux function psi, computes the neumann boundary datum g
            using bundary element method
        '''

        V = psi.function_space()
        g = Function(V)     # Neumann datum

        # compute K(psi):
        self.compute_K(psi)
        # compute Neumann boundary datum:
        g.dat.data[self.Q_dofs] = self.solve_boundary_integral_eq(psi)

        # TEST FOR g:
        g_integral = assemble(g * ds(domain=V.mesh()))   # l'integrale di g è circa 0 solo alla prima iterazione!
        print(f'Integrale di g su tutto il bordo: {g_integral}')

        return g