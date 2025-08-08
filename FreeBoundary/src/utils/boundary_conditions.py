from firedrake import *
from firedrake.mesh import plex_from_cell_list
import numpy as np
from scipy import special

def Green_function(xr,xz,yr,yz):
        '''
            Evaluates the Green function for the Grad-Shafranov operator in two given points of the poloidal plane.

            param xr: radial coordinate of point x
            param xz: height coordinate of point x
            param yr: radial coordinate of point y
            param yz: height coordinate of point y

            returns: G(x,y)
        '''
        mu0 = 4e-7 * pi
        k2 = 4*xr*yr / ((xr+yr)**2+(xz-yz)**2)
        k = sqrt(k2)
        Kk = special.ellipk(k)  # Elliptic integral of the first kind
        Ek = special.ellipe(k)  # Elliptic integral of the second kind
        return mu0 * sqrt(xr*yr) / (2*pi*k) * ( (2 - k2)*Kk - 2*Ek ) 


class JN_coupling_BCs:
    '''
        Solver for the boundary integral equation of Johnson-Nédélec.
        The solution of this equation is used as boundary condition on the artificial domain
        for the free-boundary Grad-Shafranov equation. 
    '''

    def __init__(self,params):
        '''
            params should contain:
            V - function space for the Grad-Shafranov solution "psi"
            psi - function to initialize the unknown psi, default is constant 0
            coils j - array with the toroidal current density in each coil
            coils tag - array with the mesh tag of each coil
            boundary tag - mesh tag for the artificial boundary
        '''

        V = params['V']
        self.coils_j = params['coils j']
        self.coils_tag = params['coils tag']
        self.boundary_tag = params['boundary tag']

        # Create 1D mesh on the boundary:
        m2d = V.mesh()
        self.boundary_mesh(m2d)
        self.m.init()
        
        # Function space on the boundary mesh:
        V_info = V.ufl_element()
        self.Q = FunctionSpace(self.m, V_info.family(), V_info.degree())

        # Initialize the functions psi and q
        self.psi = Function(V).interpolate(params.get('psi',Constant(0.0)))
        self.q = Function(V) 

        # Extract the needed indeces of the m2d dofs:
        self.dofs_indexes()

        # Compute functions K and V:
        self.K_func = Function(self.Q)
        self.L_func = Function(self.Q)
        self.assemble_terms()


    def boundary_mesh(self, m2d):
        '''
            "Extract" a one-dimensional mesh for the boundary of the given two-dimensional mesh.
            Nodes indexing is preserved, external facets indexing is not.
            Param m2d: Two-dimensional mesh
        '''

        # Fill "dof_coords" vector with the coordinates of the boundary nodes of the 2D mesh
        coord_func = Function(VectorFunctionSpace(m2d, "CG", 1)).interpolate(as_vector(SpatialCoordinate(m2d)))
        self.Q_dofs = FunctionSpace(m2d, "CG", 1).boundary_nodes(self.boundary_tag)
        dof_coords = coord_func.dat.data_ro[self.Q_dofs]

        # Fill "segments" with the 1D cells (segments) on the boundary
        n = len(dof_coords)
        segments = set()   # format "set" to ignore duplicates

        for i in range(n):
            # Compute the distance with the others dofs:
            dist = np.zeros(n)
            for j in range(n):
                dist[j] = np.linalg.norm(dof_coords[i]-dof_coords[j])
            dist[i] = np.inf    # set distance with itself = inf

            # Identify the indexes of the two neighbouring dofs of dof_coords[i]
            neighbour_dofs = np.argsort(dist)[:2]
            for k in neighbour_dofs:
                segments.add(tuple(sorted([i,k])))

        # Convert "segments" to array:
        segments = np.array(list(segments))

        # Build 1D boundary mesh from plex:
        plex = plex_from_cell_list(1, segments, dof_coords, comm=m2d.comm)
        self.m = Mesh(plex, dim=1, reorder = False)


    def trace(self,f):
        '''
            Extract the trace of funciton f on the boundary mesh.
            param f: function defined on the 2D mesh from which "self.m" is built in "self.boundary_mesh"
        ''' 
        # Boundary nodes of the 2D mesh:
        V = f.function_space()

        # Assign values of f on boundary nodes to its trace:
        tr_f = Function(self.Q)
        tr_f.dat.data[:] = f.dat.data_ro[self.Q_dofs]
        return tr_f  
    

    def dofs_indexes(self):
        '''
            This function is needed in the constructor.
            It extracts and stores in "V_dofs" the indexes of the dofs that lie on the boundary triangles and on the coils.
            Over these nodes the green function will be evaluated in the computation of the integral functions K and L.
            This operation is done in order to avoid computing the green function over the whole domain since it is only needed
            at the boundary and in the coils + its gradient is needed at the boundary.
        '''
        # Extract the dof of the boundary mesh cells:
        V = self.psi.function_space()
        m2d = V.mesh()
        bdry_cells = m2d.topology.exterior_facets.facet_cell_map.values    #idx of boundary cells
        V_dofs = set()
        cell_nodes = V.cell_node_map()
        for dofs in cell_nodes.values[bdry_cells]:      # nodes format of kind: [[i,j,k]]
            for dof in dofs[0]:
                V_dofs.add(dof)

        # Per qualche motivo, non tutti i dof di "Q_dofs" sono inclusi in "V_dofs".
        # Per metterci una pezza -> faccio "a mano" l'aggiunta dei dof mancanti:
        for dof in self.Q_dofs:
            V_dofs.add(dof)

        # Add to "V_dofs" the dofs in coils:
        for tag in self.coils_tag:
            coil_cells = m2d.cell_subset(tag).indices
            for dofs in cell_nodes.values[coil_cells]:
                for dof in dofs:        # Non capisco perché qua non ci vada [0] come sopra, ma così funziona
                    V_dofs.add(dof)

        self.V_dofs = np.array(list(V_dofs))


    def assemble_terms(self):
        '''
            Define the integral functions K and L. Assemble the matrix M that allows to compute q from V(q).
            The function L and the matrix M are fixed, while K needs to be updated at every iteration.
        '''

        # Extract boundary dofs indexes (Q dofs indexes)
        V = self.psi.function_space()
        m2d = V.mesh()
        x,y = SpatialCoordinate(m2d)
        mu0 = 4e-7 * pi

        # Dof coordinates for function evaluation:
        coord_func = Function(VectorFunctionSpace(m2d, "CG", 1)).interpolate(as_vector(SpatialCoordinate(m2d)))
        dof_coords = coord_func.dat.data_ro[:]

        # Matrix M for the computation of q:
        N = len(self.Q_dofs)
        self.M = np.zeros((N,N))

        self.G_list = []  # List containing the Green function for each position of the point source 

        for i in range(N):
            dof_idx = self.Q_dofs[i]    # index of the i-esim dof of Q in the 2d mesh
            X = dof_coords[dof_idx]

            # Fixed x, define G(x,y) on the boundary elements only
            G = Function(V).zero() 
            for j in self.V_dofs:
                Y = dof_coords[j]
                G.dat.data[j] = Green_function(X[0],X[1],Y[0],Y[1])

            # Fix inf value in the point source:
            G.dat.data[dof_idx] = 1e3 # Values of G are of the order of 10^-6/10^-7
            self.G_list.append(G)

            # Problema: K_func contiene tutti "nan" ed L_func tutti 0. Quindi potrebbe esserci un problema nella Green function!!
            #print("\nGreen function values for x on the boundary:")
            #print(G.dat.data_ro[self.V_dofs])
            #print(f"Valore di G nel dof della point source: {G.dat.data_ro[dof_idx]}")
            #print(f'Posizione della point source: [{X[0]},{X[1]}]')

            # Compute K value at dof i:
            n = FacetNormal(m2d)    # Outward pointing normal -> "-" needed to have inward pointing
            self.K_func.dat.data[i] = assemble( - 1 / (mu0 * x) * dot(grad(G),n) * self.psi * ds(self.boundary_tag,domain=m2d) )

            # Compute L values at dof i:
            for k in range(len(self.coils_j)):
                self.L_func.dat.data[i] += self.coils_j[k] * assemble( G * dx(self.coils_tag[k], domain=m2d))

            # Assemble matrix M:
            for l in range(N):
                dof_idx_2 = self.Q_dofs[l]
                Xl = dof_coords[dof_idx_2]

                # compute the sum of the distances between Xl and neighbouring dofs
                dist = np.zeros(N)
                for m in range(N):
                    Xm = dof_coords[self.Q_dofs[m]]
                    dist[m] = np.linalg.norm(Xl-Xm)
                dist[l] = np.inf
                sum_length = np.sum(np.sort(dist)[:2])
                
                # Fill the matrix M:
                self.M[l,i] = sum_length / 2 * G.at(Xl[0],Xl[1])


    def update_function_K(self):
        '''
            Update integral function K based on the value of the magnetic flux psi.
        '''

        V = self.psi.function_space()
        m2d = V.mesh()
        x,y = SpatialCoordinate(m2d)
        n = FacetNormal(m2d)

        mu0 = 4e-7 * pi

        # Assemble the new integral function K:
        for i in range(len(self.Q_dofs)):
            # G = Function(V).assign(self.G_list[i])
            self.K_func.dat.data[i] = assemble( - 1 / (mu0 * x) * dot(grad(self.G_list[i]),n) * self.psi * ds(self.boundary_tag,domain=m2d) )


    def solve_boundary_integral_eq(self):
        '''
            Solves the boundary integral equation for V(q), solves the linear system to compute q from V(q).
        '''

        # Extract the trace of psi:
        psi_trace = Function(self.Q)
        psi_trace.dat.data[:] = self.psi.dat.data_ro[self.Q_dofs]
        
        # Solve the problem to compute V(q):
        #v = TrialFunction(self.Q)
        #p = TestFunction(self.Q)
        #a = v * p * dx(domain=self.m)
        #L = (0.5 * psi_trace + self.K_func - self.L_func) * p * dx(domain=self.m)
        
        Vq = Function(self.Q).assign(0.5 * psi_trace + self.K_func - self.L_func)
        #solve(a==L,Vq)

        #print("Vq values before solve:") -> in Vq sono tutti nan!!
        #print(Vq.dat.data_ro[:])

        # Solve the linear system to obtain q from V(q):
        N = len(self.Q_dofs)
        #V_vec = np.zeros(N)
        #for i in range(N):
        #    dof_idx = self.Q_dofs[i]
        #    V_vec[i] = Vq.dat.data_ro[dof_idx]
        q_vec = np.linalg.solve(self.M, Vq.dat.data_ro[:])

        # Update q:      
        for i in range(N):
            dof_idx = self.Q_dofs[i]
            self.q.dat.data[dof_idx] = q_vec[i]  # q!=0 only on boundary dofs

        # TEST FOR MATRIX:
        for i in range(N):
            computed_Vq = assemble(self.G_list[i] * self.q * ds(self.boundary_tag))
            print(f'Value of V(q) = {Vq.dat.data_ro[i]}, Computed integral: {computed_Vq}')

        # TEST FOR q:
        #q_integral = assemble(self.q * ds(self.boundary_tag))   # l'integrale di q è circa 0 solo alla prima iterazione!
        #print(f'Integrale di q su tutto il bordo: {q_integral}')


    def linear_form(self,psi_old):

        self.psi.interpolate(psi_old)
        self.update_function_K()
        self.solve_boundary_integral_eq()

        V = self.psi.function_space()
        phi = TestFunction(V)
        return self.q * phi * ds(self.boundary_tag)
