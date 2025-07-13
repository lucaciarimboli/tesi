from firedrake import *

def form_a(mesh, x, phi, psi):
    """
    Define the bilinear form a in the weak formulation of the Grad-Shafranov problem.
    Parameters:
    - mesh: The mesh on which the problem is defined.
    - x: The radial coordinate (scalar).
    - phi: The test function.
    - psi: The trial function.

    Returns
    - a: the bilinear form
    """

    mu0 = 4e-7 * pi  # Permeability of free space (in IS units)
    a = ( 1 / (mu0 * x) * dot(grad(psi), grad(phi)) ) * dx(domain=mesh)
    return a


def form_b(mesh, x, G, phi, psi_N, plasma_mask, vacuum_tag):
    """
    Define the form b in the weak formulation of the Grad-Shafranov problem.
    Parameters:
    - mesh: The mesh on which the problem is defined.
    - x: The radial coordinate (scalar).
    - G(x,psi_N): lambda function representing the toroidal current density in the plasma
    - phi: The test function.
    - psi_N: The normalized poloidal magnetic flux.
    - plasma_mask: A smoothed Heavyside mask which activated in the plasma region
    - vaccum_tag: Mesh tag for the vacuum region, which contains the plasma region

    Returns
    - b: the nonlinear form
    """

    # Plasma current contribution:
    b = plasma_mask *  G(x,psi_N) * phi * dx(vacuum_tag, domain=mesh)
    return b


def form_c(mesh, phi, j_coils, coils_tags, j_cv, vessel_tag = None):
    """
    Define the linear form c in the weak formulation of the Grad-Shafranov problem.
    Parameters:
    - mesh: The mesh on which the problem is defined.
    - phi: The test function.
    - j_cv: The toroidal current density in the vessel wall
    - j_coils: Array with the toroidal current densities in each coil
    - vessel_tag: Mesh tag for the vacuum vessel wall
    - coils_tag: Array with the mesh tags for each coil

    Returns
    - c: the linear form
    """

    # Vessel wall contribution:
    #c = j_cv * phi * dx(vessel_tag, domain=mesh)
    c = 0.0 * dx(domain=mesh)
    
    # Coil contribution::
    for i in range(len(j_coils)):
        c += j_coils[i] * phi * dx(coils_tags[i], domain=mesh)  # Coils are indexed from 1 to 12 in the mesh
    
    return c


def form_d(mesh, V, phi, psi, psi0, psi_max, plasma_mask, vacuum_tag, G, dG, dpsidn):
    """
    !!!!!!!!!!!!!!!!!!!!!!!!!
    !!! INCOMPLETE METHOD !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!
    
    Define the linear form c in the weak formulation of the Grad-Shafranov problem.
    Parameters:

    Returns
    - d: the bilinear form
    """

    # Gateaux derivative of psi_N w.r.t. psi:
    psi_0 = Constant(psi0)
    psi_1 = Constant(psi_max)

    e = Function(V)
    e.interpolate( (1/(psi0 - psi_max))**2 * ((psi_0-psi_1)*psi + (psi-psi_0)*psi_1 + (psi_1-psi)*psi_0) )

    # Surface integral term:
    d = plasma_mask * dG(x,psi_N) * e * phi * dx(vacuum_tag, domain = mesh)

    # Boundary term:
    '''
    Add the integral over the plasma boundary of G(x,psi_N) * [psi(x0) - psi] * 1/dpsidn phi ...
    - Come calcolo un integrale di linea?
    - Come calcolo la derivata normale di psi rispetto ad un'isolinea?
    - Si può inserire la trial function psi come valutata in un punto? Come si fa?
    '''

    return d
