from firedrake import *

def Picard_varf(mesh, x, G, phi, psi, psi_N, plasma_mask, j_cv, I_coils, vacuum_tag, vessel_tag, coils_tags):
    """
    Define the bilinear form for the Picard iteration in the context of a magnetohydrodynamic problem.
    Parameters:
    - mesh: The mesh on which the problem is defined.
    - x: The radial coordinate (scalar).
    - G: A function representing the magnetic field or a related quantity.
    - phi: The test function.
    - psi: The trial function.
    - psi_old: The previous iteration of the trial function.
    - j_cv: The current density in the vessel.
    - I_coils: A list of current values in the coils.
    - vacuum_tag: The tag for the vacuum region in the mesh.
    - vessel_tag: The tag for the vessel region in the mesh.
    - coils_tags: A list of tags for the coil regions in the mesh.

    Returns:
    - a: The bilinear form representing the system of equations for Picard (fixed-point) iteration.
    """

    # Define the bilinear form:
    a  = (dot(grad(psi), grad(phi)) + (1/x) * Dx(psi, 0) * phi) * dx(domain=mesh) \
        - plasma_mask * G(x,psi_N) * phi * dx(vacuum_tag, domain=mesh) \
        - j_cv * phi * dx(vessel_tag, domain=mesh)
    
    # Add the contribution from the coils:
    for i in range(len(I_coils)):
        a -= I_coils[i] * phi * dx(coils_tags[i], domain=mesh)  # Coils are indexed from 1 to 12 in the mesh
    
    return a