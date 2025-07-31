import gmsh

def get_tags(geometry, params=None):
    """
    Defines in a dictionary the tags for the mesh entities needed by the solver.
    
    Parameters:
    - geometry: The type of tokamak geometry, which might be customized and built by "geometry.py" or be an existing one ("ITER").
    - params: A dictionary containing simulation parameters, such as coils and limiter points, if the geometry is custom.

    Returns:
    - A dictionary containing the tags for the mesh entities needed by the solver.
    """
    
    tags = {
        'boundary': None,
        'vacuum': None,
        'vessel': None,
        'coils':  [],
        'limiter': None,
        'limiter_pts': [],
    }
    
    if geometry == "ITER":
        # ITER geometry tags
        tags['boundary'] = 15
        #tags['vacuum'] = 13
        tags['limiter'] = 16
        tags['inside limiter'] = 2
        tags['coils'] = [3,4,5,6,7,8,9,10,11,12,13,14]
        tags['Gamma Neumann'] = 26
        tags['Gamma Dirichlet'] = 27

    elif geometry == "custom":
        tags['boundary'] = 0
        tags['vacuum'] = 2
        tags['vessel'] = 3
        if( params.get("limiter_line", None) is not None ):
            tags['limiter'] = 4
        
        tags['coils'], tags['limiter_pts'] = get_tags_from_file()

    else:
        raise ValueError(f"Geometry '{geometry}' is not recognized.")
    
    return tags


def get_tags_from_file():
    """
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

        return coils_tags, limiter_pts_tags


    except Exception as e:
        print(f"Error occurred while processing GMSH file: {e}")
        return None

    finally:
        gmsh.finalize()