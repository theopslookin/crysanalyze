"""
crysanalyze - A CLI tool for rapid preliminary analysis of powder XRD data
"""

__version__ = "0.1.0"

from .core.io import read_xy_data, write_xy_data
from .core.processing import find_peaks, subtract_background
from .core.indexing import index_peaks
from .core.lattice import calculate_lattice_parameters
from .core.texture import analyze_texture
from .core.plotting import plot_xrd, plot_peaks, plot_texture

__all__ = [
    'read_xy_data',
    'write_xy_data',
    'find_peaks',
    'subtract_background',
    'index_peaks',
    'calculate_lattice_parameters',
    'analyze_texture',
    'plot_xrd',
    'plot_peaks',
    'plot_texture'
] 