from firedrake import *
import numpy as np
from scipy import special

def Green_function(xr,xz,yr,yz):
        '''
            Evaluates the Green function for the Grad-Shafranov operator in two given points of the poloidal plane.

            param xr: radial coordinate of point x
            param xz: height coordinate of point x
            param yr: radial coordinate of point y
            param yz: height coordinate of point y

            returns: G(x,y)
        '''
        mu0 = 4e-7 * pi
        k2 = 4*xr*yr / ((xr+yr)**2+(xz-yz)**2)
        k = sqrt(k2)
        Kk = special.ellipk(k)  # Elliptic integral of the first kind
        Ek = special.ellipe(k)  # Elliptic integral of the second kind
        return mu0 * sqrt(xr*yr) / (2*pi*k) * ( (2 - k2)*Kk - 2*Ek ) 


def matrix_diagonal(len_I, ri):
    '''
        Compute the contribution of nodes neighhboring to the point source position
        for the matrix for the extraction of the Neumann data.
        Contribution are computed considering the asymptotic behavior of the Green function
        close to the point source, and considering the Neumann datum piecewise linear.

        @param len_I: boundary element length
        @param ri: radial coordinate of i-esim dof
    '''
    m = 1/(4*pi) * (len_I * (np.log(8*ri / len_I) - 1/2))
    return m

def matrix_close(len_I, ri):
    '''
        Compute the contribution of nodes neighhboring to the point source position
        for the matrix for the extraction of the Neumann data.
        Contribution are computed considering the asymptotic behavior of the Green function
        close to the point source, and considering the Neumann datum piecewise linear.

        @param len_I: boundary element length
        @param ri: radial coordinate of i-esim dof
    '''
    m = 1/(4*pi) * (len_I * (np.log(8*ri / len_I) - 3/2))
    return m


def matrix_far(len_I, Gj, rj):
    '''
        Compute the contribution of nodes far to the point source position
        for the matrix for the extraction of the Neumann data.
        Contribution are computed considering linear piecewise Neumann datum
        and Green function, i.e. using trapezoidal quadrture.

        @param len_I: boundary element length
        @param G: Green function
        @param dof_j: j-esim dof index (in the 2d mesh)
        @param rj: radial coordinate of j-esim dof

        @return m: contribution of boundary element i to j-esim dof component of matrix M
    '''
    m = 1/2 * len_I * Gj / rj
    return m


def K_neighborhood_integral(len_I, source_r, psi_source, psi_neighbor):
    '''
        Compute the line integrale contribution of the neighborhood of 
    '''

    I = len_I/(16*pi) * (psi_source - psi_neighbor) + \
        len_I/(8*pi) * np.log(8*source_r/len_I) * (psi_source + psi_neighbor)
  
    return I