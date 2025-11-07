from firedrake import *
import numpy as np

def compute_j_coils(m, tags, I):
    """
    Compute the current density flowing in each coil. The current density j is given by I / S,
    the surface is computed by Lebesgue measure using numerical integration. 

    @param mesh: The mesh on which the problem is defined.
    @param tags: The mesh tags of each coil
    @param I: Currents in Ampére in each coil

    @return j_coils: The current density for each coil [A/m^2]
    """
    
    j_coils = []
    coils_number = len(I)

    for n in range(coils_number):

        # Compute coil section size:
        S = assemble(Constant(1.0) * dx(tags[n],domain=m))

        # Append coil current density
        j_coils.append( I[n] / S)
    
    return np.array(j_coils)
        