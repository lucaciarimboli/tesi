from firedrake import *

def Picard_varf(mesh, x, G, phi, psi, psi_N, plasma_mask, j_cv, j_coils, vacuum_tag, vessel_tag, coils_tags):
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
    #a  = (dot(grad(psi), grad(phi)) + (1/x) * Dx(psi, 0) * phi) * dx(domain=mesh) \
    #    - plasma_mask *  G(x,psi_N) * phi * dx(vacuum_tag, domain=mesh) \
    #    - j_cv * phi * dx(vessel_tag, domain=mesh)

    mu0 = 4e-7 * pi  # Permeability of free space (in SI units)

    a = ( 1 / (mu0 * x) * dot(grad(psi), grad(phi)) ) * dx(domain=mesh) \
        - plasma_mask *  G(x,psi_N) * phi * dx(vacuum_tag, domain=mesh) \
        - j_cv * phi * dx(vessel_tag, domain=mesh)
    
    # Add the contribution from the coils:
    for i in range(len(j_coils)):
        a -= j_coils[i] * phi * dx(coils_tags[i], domain=mesh)  # Coils are indexed from 1 to 12 in the mesh
    
    return a



def Newton_varf(mesh, x, j, phi, psi, psi_N, psi_old, psi_denom, plasma_mask, j_cv, j_coils, vacuum_tag, vessel_tag, coils_tags):
    """
    Define the bilinear form for the Newton iteration in the context of a magnetohydrodynamic problem.
    Parameters:
    - mesh: The mesh on which the problem is defined.
    - x: The radial coordinate (scalar).
    - j: The poloidal current density in the plasma region
    - phi: The test function.
    - psi: The trial function for the flux.
    - psi_N: The normalized flux at previous iteration.
    - psi_old: The flux at the previous iteration,
    - psi_denom: psi_max - psi0. It is the denominator of the normalized psi.
    - plasma_mask: The mask that identifies the plasma region inside the vacuum chamber.
    - j_cv: The current density in the vessel walls.
    - j_coils: A list of current densities in coils.
    - vacuum_tag: The tag for the vacuum region in the mesh.
    - vessel_tag: The tag for the vessel region in the mesh.
    - coils_tags: A list of tags for the coil regions in the mesh.

    Returns:
    - a: The bilinear form representing the system of equations for Newton's method.
    """

    mu0 = 4e-7 * pi  # Permeability of free space (in SI units)

    # Define the Jacobian of the poloidal j:
    #dj_dpsi = - 1 / psi_denom * diff( j(x,psi_N), psi_N)
    '''
    r0 = 6.2
    alpha = 2.0
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    dG_dpsi = lambda_ * (beta * x / r0 + (1 - beta) * r0 / x) * gamma * alpha * psi_N**(alpha - 1) * (1 - psi_N**alpha)**(gamma - 1) / psi_denom
    '''
    '''
    # Define the bilinear form:
    a = ( 1 / (mu0 * x) * dot(grad(psi), grad(phi)) ) * dx(domain=mesh) \
        - plasma_mask * dj_dpsi * psi * phi * dx(vacuum_tag, domain=mesh)

    # Add the contribution from the r.h.s
    a -= - plasma_mask * dj_dpsi * psi_old * phi * dx(vacuum_tag, domain=mesh) \
        + plasma_mask *  j(x,psi_N) * phi * dx(vacuum_tag, domain=mesh) \
        + j_cv * phi * dx(vessel_tag, domain=mesh)
    
    # Add the contribution from the coils:
    for i in range(len(j_coils)):
        a -= j_coils[i] * phi * dx(coils_tags[i], domain=mesh)  # Coils are indexed from 1 to 12 in the mesh
    
    return a
    '''
    # Define residual as function of the previous step psi.
    F = ( 1 / (mu0 * x) * dot(grad(psi), grad(phi)) ) * dx(domain=mesh) \
        - plasma_mask *  j(x,psi_N) * phi * dx(vacuum_tag, domain=mesh) \
        - j_cv * phi * dx(vessel_tag, domain=mesh)

    for i in range(len(j_coils)):
        F -= j_coils[i] * phi * dx(coils_tags[i], domain=mesh)

    # Define the Jacobian:
    dF_dpsi = derivative(F,psi)
    dF_dpsi_N = derivative(F,psi_N)
    J = dF_dpsi - (1 / psi_denom) * dF_dpsi_N

    return J, F

    
