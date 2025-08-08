import sys
import time
sys.path.append(".")

from solver import GradShafranovSolver
from firedrake import Constant
import numpy as np


#--------------------------------------------------#
#        DEFINE SPLINES FOR CUSTOM TOKAMAK         #
#--------------------------------------------------#
# Outer boundary (1.0, 1.0) square domain by default

# Inner and outer vessel walls as circles:
def circle_points(x0, y0, r, n=100):
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return [(x0 + r * np.cos(t), y0 + r * np.sin(t)) for t in theta]

x0 = 0.7
y0 = 0.5
radius = 0.25
thick = 0.03

vessel_outer = circle_points(x0, y0, radius + thick, n=100)
vessel_inner = circle_points(x0, y0, radius, n=100)

lim_line = circle_points(x0, y0, radius - 2*thick, n=100)

# Two squared coils:
def rectangle_points(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]

coils_adapted_format = [
    [0.15, 0.25, 0.45, 0.55],
    [0.25, 0.35, 0.7, 0.8],
]

coil1 = rectangle_points(coils_adapted_format[0][0], coils_adapted_format[0][1],
                         coils_adapted_format[0][2], coils_adapted_format[0][3])
coil2 = rectangle_points(coils_adapted_format[1][0], coils_adapted_format[1][1],
                         coils_adapted_format[1][2], coils_adapted_format[1][3])

#--------------------------------------------------#
#                CURRENT PROFILE                   #
#--------------------------------------------------#

# Plasma current density profile G:
def G(R, psi_norm):
    r0 = 6.2
    alpha = 2.0
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * abs(1 - psi_norm**alpha) ** gamma

# Derivative of G w.r.t. psiN for Newton iterations:
def dGdpsiN(R, psiN):
    r0 = 6.2
    alpha = 2.0
    beta = 0.5978
    gamma = 1.395
    lambda_ = 1.365461e6
    return lambda_ * (beta * R / r0 + (1 - beta) * r0 / R) * gamma * abs(1 - psiN**alpha) ** (gamma-1) * (-alpha*psiN**(alpha-1))

#--------------------------------------------------#
#                 SET PARAMETERS                   #
#--------------------------------------------------#

params = {
    "geometry": "ITER",
    #"boundary": default,
    #"vessel_outer_wall": vessel_outer,
    #"vessel_inner_wall": vessel_inner,
    #"coils": [coil1, coil2],
    #"limiter_pts": None,
    #"limiter_line": lim_line,
    #"mesh_size_min": 0.005,
    #"mesh_size_max": 0.01,
    #"limiter_mesh_size": 0.001,
    #"dist_from_limiter": 0.1,
    #"coils_mesh_size": 0.01,
    #"dist_from_coils": 0.1,
    #"I": [-6.705e5, 1.373e4, 2.133e6, 1.432e6, -3.774e5, -6.172e5, -1.885e6, -2.359e6, -2.124e6, -1.836e6, -3.491e6, -2.04e6],
    "I": [-8.208e5, -8.621e4, 2.783e6, 1.703e6, -6.491e5, -7.971e5, -2.026e6, -2.508e6, -2.15e6, -1.874e6, -3.607e6, -2.303e6],
    "j_cv": 0,                # Vessel wall current density
    "function_space_family": "P",
    "function_space_degree": 1,
    "max_iterations": 1000,
    "tolerance": 1e-5,
    "verbose": True,
    "G": G,
    # Initial guess (can be a Constant or a Firedrake Function)
    "initial_guess": Constant(1e-4),
    "algorithm": "Picard"
}


# Build the mesh based on the geometry parameters provided
# TO BE ADJUSTED !!!!
if __name__ == "__main__":

    print("\nBuilding mesh based on custom geometry parameters...")
    path = "./meshes/custom_tokamak.msh"
            
    geometry_params = {
        "boundary": self.params.get("boundary", [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]), 
        "outer_wall": self.params.get("vessel_outer_wall", None),
        "inner_wall": self.params.get("vessel_inner_wall", None),
        "coils": self.params.get("coils", []),
        "limiter": self.params.get("limiter", None),
        "mesh_size_min": self.params.get("mesh_size_min", 0.01),
        "mesh_size_max": self.params.get("mesh_size_max", 0.05),
        "limiter_mesh_size": self.params.get("limiter_mesh_size", 0.02),
        "limiter_dist_max": self.params.get("dist_from_limiter", 0.0),
        "coils_mesh_size": self.params.get("coils_mesh_size", self.params.get("mesh_size_min")),
        "coils_dist_max": self.params.get("dist_from_coils", 0.0),
    }

    # Generate the mesh using the geometry parameters
    generate_mesh(geometry_params, path)
    self.m = Mesh(path, dim = 2, distribution_parameters={"partition": False}, reorder = True)
