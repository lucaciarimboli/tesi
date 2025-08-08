from firedrake import *

def delta_point(V,x0_idx):
    '''
       Define a Dirac delta function on the finite element space V at the point x0.
       If x0_idx contains multiple indices, than this function defines the sum of the
       dirac delta functions at those points.

       @param V: function space
       @param x0_idx: index/indeces of the nodes over which define the delta.
        indexing refers to the mesh over which V is built.
    '''

    # Define delta_x0 function: 1e6 in node x0 and 0 elsewhere
    delta_x0 = Function(V) # initialize with all zeros
    delta_x0.dat.data[x0_idx] = 1e6
    return delta_x0


def delta_line(d):
    '''
        Define a smoothed Dirac delta on a closed line.
        @param d: unsigned distance from the line function
        @TODO: adjust the smoothness so that it depends on the mesh size
    '''
    V = d.function_space()
    delta_line = Function(V).interpolate(conditional(d<2, 0.5*(1 + cos(pi*d/2)) , Constant(0.0)))
    return delta_line


def heaviside(d):
    '''
        Define a smoothed Heaviside function on a closed line.
        @param d: signed distance from the line function
        @TODO: adjust the smoothness so that it depends on the mesh size
    '''
    V = d.function_space()
    epsilon = 0.05
    heav = Function(V).interpolate(0.5 - 0.5 * tanh(d/epsilon))
    return heav