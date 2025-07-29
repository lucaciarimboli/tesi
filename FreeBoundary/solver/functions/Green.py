from firedrake import *
import numpy as np
from scipy import special

def k(xr,xz,yr,yz):
    return sqrt(4*xr*yr / ((xr+yr)**2+(xz-yz)**2))      


def Green_function_2D(V):

        m = V.mesh()
        x,y = SpatialCoordinate(m)

        


        # Let xx = [x_r,x_z] be the variable
        # Mettere le definizioni di k, K, E e tutta G in un file a parte "Green.py"
        # Le coordinate (x,y) della mesh sono la y, mentre la x sia un array di due elementi in input (xx[0],xx[1])
        # !!!!! placeholder !!!!
        return 0