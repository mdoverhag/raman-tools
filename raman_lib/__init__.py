"""
Raman Analysis Library

A library for processing and analyzing Raman spectroscopy data,
with a focus on SERS (Surface Enhanced Raman Scattering) multiplex deconvolution.
"""

from .molecules import (
    RAMAN_MOLECULES,
    CONJUGATE_COLORS,
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

from .deconvolution import (
    deconvolve_nnls
)

from .plotting import (
    plot_reference,
    plot_sample,
    plot_normalization,
    plot_deconvolution,
    plot_deconvolution_original_scale,
    plot_deconvolution_boxplots,
    plot_peak_intensity_histogram
)

from .workflow import (
    build_reference_dict,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples,
    extract_peak_intensities,
    print_experiment_summary
)

__version__ = "0.1.0"

__all__ = [
    # Molecules
    "RAMAN_MOLECULES",
    "CONJUGATE_COLORS",
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

    # Deconvolution
    "deconvolve_nnls",

    # Plotting
    "plot_reference",
    "plot_sample",
    "plot_normalization",
    "plot_deconvolution",
    "plot_deconvolution_original_scale",
    "plot_deconvolution_boxplots",
    "plot_peak_intensity_histogram",

    # Workflow
    "build_reference_dict",
    "load_and_process_reference",
    "load_and_process_sample",
    "normalize_and_deconvolve_samples",
    "extract_peak_intensities",
    "print_experiment_summary",
]
