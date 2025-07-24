from firedrake import *
from firedrake.mesh import plex_from_cell_list
import numpy as np

class JN_BCs:
    '''
        Solver for the boundary integral equation of Johnson-Nédélec.
        The solution of this equation is used as boundary condition on the artificial domain
        for the free-boundary Grad-Shafranov equation. 
    '''

    def __init__(self,params):

        # Create 1D mesh on the boundary:
        self.boundary_mesh(params['mesh'])
        self.m.init()


    def boundary_mesh(self, m2d):
        '''
            "Extract" a one-dimensional mesh for the boundary of the given two-dimensional mesh.
            Nodes indexing is preserved, external facets indexing is not.
            Param m2d: Two-dimensional mesh
        '''

        # Fill "dof_coords" vector with the coordinates of the boundary nodes of the 2D mesh
        coord_func = Function(VectorFunctionSpace(m2d, "CG", 1)).interpolate(as_vector(SpatialCoordinate(mesh)))
        dofs = DirichletBC(FunctionSpace(m2d, "CG", 1), 0.0, "on_boundary").nodes
        dof_coords = coord_func.dat.data_ro[dofs]

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
            # Convert "segments" in array:
            segments = np.array(list(segments))

        # Build 1D boundary mesh from plex:
        plex = plex_from_cell_list(1, segments, dof_coords, comm=m2d.comm)
        self.m = Mesh(plex, dim=1, reorder = False)

    def trace(f):
        '''
            
        ''' 
