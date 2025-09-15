"""
AbacusSummit Intrinsic Alignment Analysis Package

This package provides tools for generating LRG catalogs from AbacusSummit
simulation data using HOD models, including intrinsic alignment measurements
through tidal shear calculations.
"""

__version__ = "1.0.0"
__author__ = "LRG HOD Pipeline"

from .lrg_hod_pipeline import LRGHODPipeline
from .efficient_tidal_shear import compute_efficient_tidal_shear

__all__ = [
    "LRGHODPipeline",
    "compute_efficient_tidal_shear",
]