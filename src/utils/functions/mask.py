from firedrake import *
import numpy as np

def delta_point(V,x0_idx):
    '''
    Define a Dirac delta function on the finite element space V at the point x0.
    If x0_idx contains multiple indices, than this function defines the sum of the
    dirac delta functions at those points.

    @param V: function space
    @param x0_idx: index/indeces of the nodes over which define the delta.
    '''

    # Define delta_x0 function: 1e6 in node x0 and 0 elsewhere
    delta_x0 = Function(V) # initialize with all zeros
    delta_x0.dat.data[x0_idx] = 1e6
    return delta_x0


def delta_FD(x,h,m=1):
    '''
    Returns the dirac delta approximation obtained as the derivative
    of the Fermi-Dirac function.
    x is the distance, h is the mesh size and m is a scaling coefficient for h.
    '''
    eps = h*m
    return 1 / eps * exp(-x/eps)/(1+exp(-x/eps))**2


def delta_line(d,h):
    '''
    Define a smoothed Dirac delta on a closed line.
    @param d: unsigned distance from the line function
    '''
    m = 1

    V = d.function_space()
    delta_line = Function(V).interpolate(conditional(d<6*h, delta_FD(d,h,m) , Constant(0.0)))
    return delta_line


def heaviside(d,h):
    '''
    Define a smoothed Heaviside function on a closed line.
    @param d: signed distance from the line function
    '''
    V = d.function_space()
    heav = Function(V).interpolate(0.5 - 0.5 * tanh(d/h))
    return heav