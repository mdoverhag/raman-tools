"""
Raman Analysis Library

A library for processing and analyzing Raman spectroscopy data,
with a focus on SERS (Surface Enhanced Raman Scattering) multiplex deconvolution.
"""

from .molecules import (
    RAMAN_MOLECULES,
    get_peak,
    get_color,
    get_all_molecules
)

from .io import (
    load_spectrum,
    load_spectra,
    create_output_dir,
    ensure_output_subdir
)

__version__ = "0.1.0"

__all__ = [
    # Molecules
    "RAMAN_MOLECULES",
    "get_peak",
    "get_color",
    "get_all_molecules",

    # I/O
    "load_spectrum",
    "load_spectra",
    "create_output_dir",
    "ensure_output_subdir",
]
