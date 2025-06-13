from firedrake import *

def compute_j_coils(mesh, tags, I):
    
    j_coils = []

    coils_number = len(I)

    for n in range(0,coils_number - 1):

        # Compute coil section size:
        S = assemble(Constant(1.0) * dx(tags[n],domain=mesh))

        # Append coil current density
        j_coils.append( I[n] / S)
    
    return j_coils
        