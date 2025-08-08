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

    def __init__(self,V,limiter_tag):
        '''
            @brief Constructor.

            @param V: Function Space
            @param limiter_tag: Mesh tag for the limiter nodes
            @TODO modify so that the constructor handles whatever initial condition
            (means remove initalization here and take care of it in the solver)
        '''

        # Initialise maks for the plasma region and its boundary
        self.domain_mask = Function(V).interpolate(Constant(1.0))
        self.boundary_mask = Function(V).interpolate(Constant(0.0))

        # Build a neighbors map:
        self.build_neighbors_map(V)

        # Store the list of limiter nodes:
        self.limiter = DirichletBC(V, 0.0, limiter_tag).nodes

        # Initialize all plasma values:
        self.psi0 = 0.0
        self.psi_ax = -1.0
        self.x0_idx = np.array([0])
        self.x1_idx = np.array([0])
        self.n = Function(VectorFunctionSpace(V.mesh(), "Lagrange", 1, dim=2)).interpolate(as_vector((0.0, 0.0)))


    def build_neighbors_map(self, V):
        '''
        @brief Map each dof to its neighboring dofs  
        @param V: function space (expected P1 elements)

        Builds a map that provides for each dof an array containing its
        neighboring dofs (i.e. directly connected by mesh edges) ordered counter-clockwise.
        The result is a list of np.arrays stored as a class member named "neighbors map".
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

        for idx, neighbors_set in enumerate(neighbors_map):
            neighbors_array = np.array(list(neighbors_set))

            node_coords = coord_func.dat.data_ro[idx]
            neighbor_coords = coord_func.dat.data_ro[neighbors_array]

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
        @param neighbors_map: list that associates to a node an array with the indexes of the neighboring nodes
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

    
    def identify_psi0(self,psi_data):
        '''
        @param psi: flux function

        Identifies the value at the plasma-boundary contour line as the greatest among:
            - the value of psi at every limiter point
            - the value of psi at every saddle point
        The array self.saddle_points is filled with all the saddle points that lie on the 
        plasma boundary (contour line psi = psi0), or leaved empty if the plasma boundary
        is tangent to the limiter
        '''

        limiter_data = psi_data[self.limiter]
        self.x0_idx = self.limiter[np.argmax(limiter_data)]
        self.psi0 = max(limiter_data)

        # Extract nodes s.t. psi>psi0
        internal_dof_idx = np.where(psi_data > self.psi0)[0] 

        # Identify saddle points inside the region psi>psi0
        saddle_pts = np.array([idx for idx in internal_dof_idx if self.is_saddle_point(psi_data, idx)])

        # Update psi0 if there is an X-points configuration:
        if len(saddle_pts) > 0:
            self.psi0 = max(psi_data[saddle_pts])
            self.saddle_points = saddle_pts[psi_data[saddle_pts]==self.psi0]
        else:
            self.saddle_points = np.array([])
        

    def compute_intersections(self,psi_data,dof_coords,idx):
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
        queue = {idx}
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
                elif neighbor in self.saddle_points:
                    # because "visited" at the end of the loop is the set of dofs inside inside the contour line
                    # so even if the saddle point it is not actually visited, it is added to include in the set
                    # of dofs inside the plasma boundary
                    visited.add(neighbor)
                    continue
                else:
                    queue.add(neighbor)

        return np.array(intersections), np.array(list(visited))


    def update(self,psi):
        '''
            Update the plasma boundary masks provided the current flux solution psi
            @params psi: flux function
            
            @TODO: Modify the criterion of the masks so that the width of the
                plasma boundary mask and the smoothness coeff of the plasma mask depend on the
                mesh size h.
                Currently the width is fixed to "2" and the smoothness is "epsilon"
        '''

        psi_data = psi.dat.data_ro[:]
        V = psi.function_space()
        m = V.mesh()

        # Identify the value of psi at the plasma boundary and at the magnetic axis
        self.identify_psi0(psi_data)
        self.x1_idx = np.argmax(psi_data) # magnetic axis index(es)
        self.psi_ax = max(psi_data)

        # Extract coordinates of DOFs
        coord_func = Function(VectorFunctionSpace(m, "CG", 1)).interpolate(as_vector(SpatialCoordinate(m)))
        dof_coords = coord_func.dat.data_ro[:]

        # Define an unsigned distance function from the plasma boundary:
        level_pts, inside_dofs = self.compute_intersections(psi_data, dof_coords, self.x1_idx)
        tree = cKDTree(level_pts)
        d = Function(V)
        for i, pt in enumerate(dof_coords):
            dist, _ = tree.query(pt)
            d.dat.data[i] = dist

        # Using the unsigned distance, define a mask for the plasma boundary:
        self.boundary_mask.assign(delta_line(d))

        # Using the unsigned distance, define the normal vector on the plasma boundary (defined everywhere)
        self.n = Function(VectorFunctionSpace(m, "Lagrange", 1, dim=2)).interpolate(grad(d))

        # Convert the distance function to signed, by changing sign at the dofs marked as inside the plasma:
        d.dat.data[inside_dofs] *= -1
        self.domain_mask.assign(heaviside(d))