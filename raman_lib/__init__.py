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

from .baseline import (
    als_baseline,
    apply_baseline_correction
)

from .averaging import (
    calculate_average
)

from .normalization import (
    normalize_l2,
    normalize_spectra_l2
)

from .plotting import (
    plot_reference,
    plot_sample,
    plot_normalization,
    plot_deconvolution
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

    # Baseline correction
    "als_baseline",
    "apply_baseline_correction",

    # Averaging
    "calculate_average",

    # Normalization
    "normalize_l2",
    "normalize_spectra_l2",

    # Plotting
    "plot_reference",
    "plot_sample",
    "plot_normalization",
    "plot_deconvolution",
]
