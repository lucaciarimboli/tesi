from firedrake import *

def compute_j_coils(mesh, tags, I):
    """
    Compute the current density flowing in each coil. The current density j is given by I / S,
    the surface is computed by Lebesgue measure using numerical integration. 

    Parameters:
    - mesh: The mesh on which the problem is defined.
    - tags: The mesh tags of each coil
    - I: Currents in Ampére in each coil

    Returns:
    - j_coils: The current density for each coil [A/m^2]
    """
    
    j_coils = []

    coils_number = len(I)

    for n in range(0,coils_number - 1):

        # Compute coil section size:
        S = assemble(Constant(1.0) * dx(tags[n],domain=mesh))

        # Append coil current density
        j_coils.append( I[n] / S)
    
    return j_coils
        