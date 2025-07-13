
from firedrake import *

def farfield_form(mesh,phi,psi,boundary_id,radius):
    mu0 = 4e-7 * pi  # Permeability of free space (in IS units)
    x,y = SpatialCoordinate(mesh)

    N = build_N(x,y,radius)
    single_integral = (1/mu0) * psi * N * phi * ds(boundary_id)

    return single_integral #+ double_integral

def build_N(x,y,radius):
    delta_plus = sqrt( x**2 + (radius + y)**2)
    delta_minus = sqrt( x**2 + (radius - y)**2)
    return 1/x * (1/delta_plus + 1/delta_minus - 1/radius)

'''
    In firedrake NON è possibile definire un doppio integrale come richiesto.
    Nel forum github di firedrake, un tizio ha già chiesto come poter fare e, verosimilmente, stava implementando
    il mio stesso problema: https://github.com/firedrakeproject/firedrake/discussions/3434.

    Quello che un utente suggerisce è di imporre non-local boundary conditions invece di usare BEM. 
    Potrei, a partire dal 
'''