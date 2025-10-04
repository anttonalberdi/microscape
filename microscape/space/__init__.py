# microscape/space/__init__.py
from .graph import build_edge_index
from .diffusion import explicit_diffusion_step, stable_dt_upper_bound
from .coarsen import coarsen

__all__ = [
    "Edge",
    "build_laplacians",
    "explicit_diffusion_step",
    "stable_dt_upper_bound",
    "coarsen",
]
