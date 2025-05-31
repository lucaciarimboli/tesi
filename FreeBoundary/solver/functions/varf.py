from firedrake import *

def GS_varf_Picard(x, G, phi, psi, psi_old, psi_mask, vessel_mask, coils_mask):
    # Define the bilinear form:
    a  = (dot(grad(psi), grad(phi)) + (1/x) * Dx(psi, 0) * phi) * dx \
        - psi_mask * G(x,psi_old) * phi * dx \
        - coils_mask * phi * dx \
        - vessel_mask * phi * dx
    
    # Note: Current in vessel walls and coils is already accounted for in the masks.
    return a