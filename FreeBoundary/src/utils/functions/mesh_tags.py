import gmsh
import numpy as np

def get_tags(geometry):
    """
    Defines in a dictionary the tags for the mesh entities needed by the solver.
    
    @param geometry: The type of tokamak geometry, which might be customized and built by "geometry.py" or be an existing one ("ITER").
    @param params: A dictionary containing simulation parameters, such as coils and limiter points, if the geometry is custom.

    @return A dictionary containing the tags for the mesh entities needed by the solver.
    """
    
    if geometry == "ITER":

        # This mesh does not account for current density on the vessel.
        # The tag is set as 1 (all domain except vacuum region and coils) as a placeholder.

        tags = {
            'boundary': 15,
            'vacuum': 2,
            'vessel': 1,
            'coils': np.array([3,4,5,6,7,8,9,10,11,12,13,14]),
            'limiter': 16
        }   

    elif geometry == "custom":

        coils_tags, _ = coils_tags_from_file()

        tags = {
            'boundary': 0,
            'vacuum': 1,
            'vessel': 3,
            'limiter': 16,
            'coils': coils_tags
        }

    else:
        raise ValueError(f"Geometry '{geometry}' is not recognized.")
    
    return tags


def get_tags_for_plot(geometry):
    
    if geometry == "ITER":

        # This mesh does not account for current density on the vessel.
        # The tag is set as 1 (all domain except vacuum region and coils) as a placeholder.

        coils_bdry_tags = [18, 19, 20, 21, 22, 23]

        tags = {
            'limiter': 16,
            'solenoid': 17,
            'coils': np.array(coils_bdry_tags)
        }

    elif geometry == "custom":

        _, coils_tags = coils_tags_from_file()

        tags = {
            'boundary': 0,
            'vacuum': 1,
            'vessel': 3,
            'limiter': 16,
            'coils': coils_tags
        }

    else:
        raise ValueError(f"Geometry '{geometry}' is not recognized.")
    
    return tags 


def coils_tags_from_file():
    """
    @TODO check and adjust this method functioning.
    Fai in modo che restituisca due array numpy: uno con i nodi dentro ogni coils ed uno
    con i nodi sul bordo delle coils.

    This function retrieves the tags for the mesh entities from a GMSH file.
    This is necessary if the number of coils or of limiter points is not known beforehand.

    Parameters:
    - geometry: The type of tokamak geometry, which might customized and built by "geometry.py" or be an existing one ("ITER").

    Returns:
    - A dictionary containing the tags for the mesh entities needed by the solver.
    """

    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    try:
        gmsh.open("./meshes/custom_tokamak.msh")

        physical_groups = gmsh.model.getPhysicalGroups()

        coils_tags = []
        limiter_pts_tags = []

        for dim, tag in physical_groups:

            # Look at the coils:
            if dim == 2:
                name = gmsh.model.getPhysicalName(dim, tag)
                if name.startswith("Coil_"):
                    coils_tags.append(tag)

            # Look at the limiter points if any:
            if dim == 0:
                name = gmsh.model.getPhysicalName(dim, tag)
                if name.startswith("LimiterPoint_"):
                    limiter_pts_tags.append(tag)

            # AGGIUNGI L'ESTRAZIONE DEL TAG DELLA REGIONE INTERNA AL LIMITER
            # NEL CASO DI LIMITER_LINE

        if not coils_tags:
            print("No coil tags found in the GMSH file.")

        return np.array(coils_tags)

    except Exception as e:
        print(f"Error occurred while processing GMSH file: {e}")
        return None

    finally:
        gmsh.finalize()