from firedrake import *
import numpy as np
from scipy.spatial import cKDTree
from src.utils.functions.mask import delta_line, heaviside

def interpolate_edge(p1, p2, v1, v2, c):
    t = (c - v1) / (v2 - v1)
    return (1 - t) * p1 + t * p2

class Plasma:
    '''
    Class to extract the plasma boundary from a flux function in a 2D tokamak.
    It looks for saddle points, identifies the plasma boundary and stores
    a mask for the plasma region and one for the plasma boundary
    '''

    def __init__(self,V,h,limiter_tag,vacuum_tag,plasma_domain=Constant(1.0), plasma_boundary=Constant(0.0)):
        '''
        Constructor.

        @oaram V: function space
        @param h: mesh size in the vacuum region
        @param limiter_tag: tag for the limiter boundary
        @param vacuum_tag: tag for the vacuum region
        @param plasma_domain: initial mask for the plasma region
        @param plasma_boundary: initial mask for the plasma boundary
        '''

        # Define mesh size in the vacuum region:
        self.h = h

        # Initialise maks for the plasma region and its boundary
        self.domain_mask = Function(V).interpolate(plasma_domain)
        self.boundary_mask = Function(V).interpolate(plasma_boundary)

        # Build a neighbors map:
        self.build_neighbors_map(V)

        # Store the list of limiter nodes:
        self.limiter = DirichletBC(V, 0.0, limiter_tag).nodes

        # Extract the set of dof that lie in the vacuum region:
        v = TestFunction(V)
        vacuum_vector = assemble(v * dx(vacuum_tag))
        vacuum_vector.dat.data[self.limiter] = 0.0  # exclude boundary points
        self.vacuum = np.where(vacuum_vector.dat.data_ro > 0.0)[0]

        # Initialise the set of dof that lie inside the plasma region:
        self.inside_dofs = np.array([], dtype=int)

        # Initialize all plasma values:
        self.psi0 = 0.0
        self.psi_ma = -1.0
        self.x0_idx = self.limiter[len(self.limiter) // 5]  # Take some index in limiter region
        self.x1_idx = self.vacuum[len(self.vacuum) // 5]    # Take some index in vacuum region
        
        #---------------- COORDINATES OF THE MAGNETIC AXIS ----------------#
        ''''
        coord_func = Function(VectorFunctionSpace(V.mesh(), "CG", 1)).interpolate(as_vector(SpatialCoordinate(V.mesh())))
        self.x1 = coord_func.dat.data_ro[self.x1_idx]
        print(f"Initial magnetic axis guess at index {self.x1_idx} with coordinates {self.x1}")
        print(f"Initial X-point guess at index {self.x0_idx} with coordinates {coord_func.dat.data_ro[self.x0_idx]}")
        '''
        #------------------------------------------------------------------#

        self.n = Function(VectorFunctionSpace(V.mesh(), "Lagrange", 1, dim=2)).interpolate(as_vector((0.0, 0.0)))

        # For plotting saddle points for debugging:
        self.psi_X_point = 0.0


    def build_neighbors_map(self, V):
        '''
        Builds a map that provides for each dof an array containing its
        neighboring dofs (i.e. directly connected by mesh edges) ordered counter-clockwise.
        The result is a list of np.arrays stored as a class member named "neighbors map".

        @param V: function space
        '''
        n_dofs = V.dof_count
        neighbors_map = [set() for _ in range(n_dofs)]  # list of set

        cell_dofs = V.cell_node_map().values
        for cell in cell_dofs:
            for i in range(len(cell)):
                dof1 = cell[i]
                for j in range(i+1,len(cell)):
                    dof2 = cell[j]
                    neighbors_map[dof1].add(dof2)
                    neighbors_map[dof2].add(dof1)

        # Order (counter-clockwise) neighbors order:
        self.neighbors_map = []
        m = V.mesh()
        coord_func = Function(VectorFunctionSpace(m, "CG", 1)).interpolate(as_vector(SpatialCoordinate(m)))
        self.coords = coord_func.dat.data_ro[:]

        for idx, neighbors_set in enumerate(neighbors_map):
            neighbors_array = np.array(list(neighbors_set))

            node_coords = self.coords[idx]
            neighbor_coords = self.coords[neighbors_array]

            dx = neighbor_coords[:, 0] - node_coords[0]
            dy = neighbor_coords[:, 1] - node_coords[1]
            angles = np.arctan2(dy, dx)

            sorted_indices = np.argsort(angles)
            self.neighbors_map.append(neighbors_array[sorted_indices])

    
    def is_saddle_point(self, psi_data, idx):
        '''
        Check if the node indicated with idx is a saddle point for the function f

        @param psi_data: array with function psi value on every dof
        @param idx: index of the mesh node to be checked.

        @return: True if idx-th dof is a saddle point, False otherwise
        '''
        neighbors_idx = self.neighbors_map[idx] 
        if len(neighbors_idx) < 4:
            return False

        # If grad changes sign at least 4 times, f has a saddle point in dof[idx]:
        f_node = psi_data[idx]
        f_neighbors = psi_data[neighbors_idx]
        grad_signs = np.sign(f_neighbors - f_node)

        # Eliminate zeros to count the number of sign changes:
        for i in range(len(grad_signs)):
            if grad_signs[i] == 0:
                grad_signs[i] = grad_signs[i-1] if i > 0 else grad_signs[-1]

        counter = 0
        for i in range(len(grad_signs)):
            next_i = (i + 1) % len(grad_signs)
            if grad_signs[i] * grad_signs[next_i] < 0:
                counter += 1

        return counter > 3  


    def identify_boundary(self,psi_data):
        '''
        Identifies the value at the plasma-boundary contour line as the greatest among:
            - the value of psi at every saddle point
            - the value of psi at the limiter points which are not in the shadow of a saddle

        To roughly select the limiter points which are not in the shadow of an x-point,
        the domain is divided in four quadrants w.r.t. the position of the magnetic axis.
        The limiter points that are in any quadrant that contains an x-point are considered
        in its shadow.

        @param psi_data: value of the flux function at every dof
        @TODO: More accurate criterion for the limiter points selection.
        '''

        # Look for saddle points in the vacuum region:
        saddle_points = np.array([idx for idx in self.vacuum if self.is_saddle_point(psi_data, idx)])

        if len(saddle_points) > 0:
            # Candidate X-point is the saddle with largest flux:
            self.x0_idx = saddle_points[np.argmax(psi_data[saddle_points])]
            self.psi0 = psi_data[self.x0_idx]

            #-------- FOR DEBUG --------#
            self.psi_X_point = self.psi0
            self.X_point_idx = self.x0_idx
            #---------------------------#

            # Identify limiter points NOT in the shadow of an x-point:
            lim_not_shadow = self.limiter.copy()
            x1 = self.coords[self.x1_idx]    # magnetic axis coordinates
            for sp in saddle_points:
                X = self.coords[sp]  # x-point coordinates
                quadrant = np.sign(X - x1)
                lim_not_shadow = np.array([
                    idx for idx in lim_not_shadow
                    if np.any(np.sign(self.coords[idx] - x1) != quadrant)
                ])

            if len(lim_not_shadow) > 0:
                limiter_idx = lim_not_shadow[np.argmax(psi_data[lim_not_shadow])]
                limiter_psi = psi_data[limiter_idx]

                # Compare flux value to find psi0:            
                if limiter_psi > self.psi0:
                    self.x0_idx = limiter_idx
                    self.psi0 = limiter_psi

        else:
            # If there are no saddles, take the limiter point where psi is maximum
            self.x0_idx = self.limiter[np.argmax(psi_data[self.limiter])]
            self.psi0 = psi_data[self.x0_idx]
        
            #-------- FOR DEBUG --------#
            self.psi_X_point = [self.psi0]
            #---------------------------#


    def compute_intersections(self,psi_data,dof_coords):
        '''
        Computes the points in which the contour line psi = self.psi0 intersect
        triangles edges.

        @param psi_data: array with function psi value on every dof
        @param dof_coords: array with the spatial coordinates of every dof 
        @param idx: index of the dof where the magnetic axis is located

        @return intersections: array with the coordinates of intersection points
        @return visited: array with the indexes of the dofs internal to the plasma
        '''

        # Explore mesh starting from the magnetic axis
        intersections = []
        queue = {self.x1_idx}
        visited = set()

        while queue:
            dof = queue.pop()
            visited.add(dof)    

            p1 = dof_coords[dof]
            f1 = psi_data[dof]
            c = self.psi0

            for neighbor in self.neighbors_map[dof]:

                if neighbor in visited:
                    continue

                f2 = psi_data[neighbor]
                if (f1 - c) * (f2 - c) <= 0:  # sign changes -> intersection
                    p2 = dof_coords[neighbor]
                    pt = interpolate_edge(p1, p2, f1, f2, c)
                    intersections.append(pt)
                elif neighbor == self.x0_idx:
                    # because "visited" at the end of the loop is the set of dofs inside inside the contour line
                    # so even if the saddle point it is not actually visited, it is added to include in the set
                    # of dofs inside the plasma boundary
                    visited.add(neighbor)
                    continue
                else:
                    queue.add(neighbor)

        return np.array(intersections), np.array(list(visited))
    

    def identify_magnetic_axis(self, psi_data, idx):
        '''
        Function to identify the magnetic axis, using a recursive
        hill-climbing algorithm

        @param psi_data: array with function psi value on every dof
        @param idx: index of the candidate magnetic axis node
        '''
        # psi value at candidate magnetic axis
        psi_ma = psi_data[idx]

        # max psi value among neighbors
        neighbors_idx = np.array([
            ni for ni in self.neighbors_map[idx]
            if ni not in self.limiter   # avoid crossing
        ])
        next_idx = neighbors_idx[np.argmax(psi_data[neighbors_idx])]

        # Stop search if idx is local maximum
        if psi_data[next_idx] > psi_ma:
            self.identify_magnetic_axis(psi_data, next_idx)
        else:
            self.x1_idx = idx
            self.psi_ma = psi_data[idx]


    def update(self,psi,imposed_x0_idx=None):
        '''
        Update the plasma boundary masks provided the current flux solution psi

        @params psi: flux function
        @param imposed_x0_idx: if provided, use this index as x-point location

        @return signed distance function to allow plotting the plasma boundary
        '''

        psi_data = psi.dat.data_ro[:]
        V = psi.function_space()
        m = V.mesh()

        # Identify the value of psi at the plasma boundary and at the magnetic axis
        # self.x1_idx = self.vacuum[np.argmax(psi_data[self.vacuum])]
        #self.psi_ma = psi_data[self.x1_idx]
        self.identify_magnetic_axis(psi_data, self.x1_idx)
        self.identify_boundary(psi_data)
        if imposed_x0_idx is not None:
            self.x0_idx = imposed_x0_idx
            self.psi0 = psi_data[self.x0_idx]

        # Extract coordinates of DOFs
        coord_func = Function(VectorFunctionSpace(m, "CG", 1)).interpolate(as_vector(SpatialCoordinate(m)))
        dof_coords = coord_func.dat.data_ro[:]

        # Define an unsigned distance function from the plasma boundary:
        level_pts, self.inside_dofs = self.compute_intersections(psi_data, dof_coords)
        tree = cKDTree(level_pts)
        self.d = Function(V)
        for i, pt in enumerate(dof_coords):
            dist, _ = tree.query(pt)
            self.d.dat.data[i] = dist

        # Using the unsigned distance, define a mask for the plasma boundary:
        self.boundary_mask.assign(delta_line(self.d,self.h))

        # Using the unsigned distance, define the normal vector on the plasma boundary (defined everywhere)
        #self.n = Function(VectorFunctionSpace(m, "Lagrange", 1, dim=2)).interpolate(grad(self.d))

        # Convert the distance function to signed, by changing sign at the dofs marked as inside the plasma:
        self.d.dat.data[self.inside_dofs] *= -1
        self.n = Function(VectorFunctionSpace(m, "Lagrange", 1, dim=2)).interpolate(grad(self.d))
        self.domain_mask.assign(heaviside(self.d,self.h))