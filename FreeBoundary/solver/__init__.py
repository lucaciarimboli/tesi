# File: /freeboundary/freeboundary/freeboundary/__init__.py

from .solver import GradShafranovSolver
from .functions.geometry import generate_mesh
from .functions.mesh_tags import get_tags
from .functions.varf import Picard_varf, Newton_varf
from .functions.coils import compute_j_coils

__all__ = [
    "GradShafranovSolver",
    "generate_mesh",
    "get_tags",
    "Picard_varf",
    "Newton_varf"
    "compute_j_coils"
]