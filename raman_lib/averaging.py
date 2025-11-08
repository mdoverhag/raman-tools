"""
Averaging and statistics for Raman spectroscopy data.

Calculate mean and standard deviation across multiple replicate spectra.
"""

import numpy as np


def calculate_average(spectra: list[dict]) -> dict:
    """
    Calculate average and standard deviation across multiple spectra.

    All spectra must have the same wavenumber values. This function assumes
    baseline correction has already been applied to each spectrum.

    Args:
        spectra: List of spectrum dictionaries, each with:
                 - 'wavenumbers': list of wavenumber values
                 - 'intensities': list of raw intensity values
                 - 'corrected': list of baseline-corrected intensity values
                 - 'baseline': list of baseline values

    Returns:
        Dictionary with averaged data:
        {
            'wavenumbers': [...],
            'raw_avg': [...],        # average of raw intensities
            'raw_std': [...],        # std dev of raw intensities
            'corrected_avg': [...],  # average of corrected intensities
            'corrected_std': [...],  # std dev of corrected intensities
            'baseline_avg': [...],   # average of baselines
            'count': int             # number of spectra averaged
        }

    Raises:
        ValueError: If spectra list is empty or spectra have different wavenumbers
    """
    if not spectra:
        raise ValueError("Cannot average empty list of spectra")

    # Verify all spectra have same wavenumbers
    reference_wavenumbers = spectra[0]['wavenumbers']
    n_points = len(reference_wavenumbers)

    for i, spectrum in enumerate(spectra):
        if len(spectrum['wavenumbers']) != n_points:
            raise ValueError(
                f"Spectrum {i} has {len(spectrum['wavenumbers'])} points, "
                f"expected {n_points}"
            )

    # Collect all intensities, corrected values, and baselines
    all_raw = np.array([s['intensities'] for s in spectra])
    all_corrected = np.array([s['corrected'] for s in spectra])
    all_baselines = np.array([s['baseline'] for s in spectra])

    # Calculate averages and standard deviations along axis 0 (across spectra)
    raw_avg = np.mean(all_raw, axis=0)
    raw_std = np.std(all_raw, axis=0, ddof=1)  # Use sample std dev (N-1)

    corrected_avg = np.mean(all_corrected, axis=0)
    corrected_std = np.std(all_corrected, axis=0, ddof=1)

    baseline_avg = np.mean(all_baselines, axis=0)

    return {
        'wavenumbers': reference_wavenumbers,
        'raw_avg': raw_avg.tolist(),
        'raw_std': raw_std.tolist(),
        'corrected_avg': corrected_avg.tolist(),
        'corrected_std': corrected_std.tolist(),
        'baseline_avg': baseline_avg.tolist(),
        'count': len(spectra)
    }
