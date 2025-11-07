import gmsh
import numpy as np

def get_tags(geometry):
    """
    Defines in a dictionary the tags for the mesh entities needed by the solver.
    
    @param geometry: The type of tokamak geometry.
    @return A dictionary containing the tags for the mesh entities needed by the solver.
    """
    
    if geometry == "ITER":

        # This mesh does not account for current density on the vessel.
        # The tag is set as 1 (all domain except vacuum region and coils) as a placeholder.

        tags = {
            'boundary': 15,
            'dirichlet_boundary': 27,
            'neumann_boundary': 26,
            'vacuum': 2,
            'vessel': 1,
            'coils': np.array([3,4,5,6,7,8,9,10,11,12,13,14]),
            'limiter': 16
        }  

    elif geometry == "CompassU":

        # This mesh does not account for current density on the vessel.
        # The tag is set as 1 (all domain except vacuum region and coils) as a placeholder.

        tags = {
            'boundary': 18,
            'dirichlet_boundary': 20,
            'neumann_boundary': 19,
            'vacuum': 2,
            'vessel': 1,    
            'coils': np.array([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]),
            'limiter': 21
        } 

    else:
        raise ValueError(f"Geometry '{geometry}' is not recognized.")
    
    return tags


def get_tags_for_plot(geometry):
    """
    Defines in a dictionary the tags for the mesh entities needed for plotting.

    @param geometry: The type of tokamak geometry.
    @return A dictionary containing the tags for the mesh entities needed for plotting.
    """
    
    if geometry == "ITER":

        # This mesh does not account for current density on the vessel.
        # The tag is set as 1 (all domain except vacuum region and coils) as a placeholder.

        coils_bdry_tags = [18, 19, 20, 21, 22, 23]

        tags = {
            'limiter': 16,
            'solenoid': 17,
            'coils': np.array(coils_bdry_tags)
        }

    elif geometry == "CompassU":

        # This mesh does not account for current density on the vessel.
        # The tag is set as 1 (all domain except vacuum region and coils) as a placeholder.

        coils_bdry_tags = [22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]

        tags = {
            'limiter': 21,
            'coils': np.array(coils_bdry_tags)
        }

    else:
        raise ValueError(f"Geometry '{geometry}' is not recognized.")
    
    return tags 