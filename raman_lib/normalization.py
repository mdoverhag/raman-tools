"""
Normalization methods for Raman spectroscopy data.

Implements L2 normalization within a specified wavenumber range.
"""

import numpy as np


def normalize_l2(
    spectrum: dict,
    wavenumber_range: tuple[float, float]
) -> dict:
    """
    Apply L2 normalization to a spectrum within a specified wavenumber range.

    L2 normalization divides the spectrum by its Euclidean norm (L2 norm)
    calculated over the specified wavenumber range. This makes spectra
    comparable regardless of overall intensity differences.

    Args:
        spectrum: Dictionary with 'wavenumbers' and 'corrected_avg' (or similar intensity key)
        wavenumber_range: Tuple of (min_wavenumber, max_wavenumber) for normalization

    Returns:
        Dictionary with:
        - 'wavenumbers': same as input
        - 'normalized': L2-normalized intensities
        - 'original': original intensities (for reference)
        - 'norm_factor': the L2 norm value used for normalization
        - 'norm_range': the wavenumber range used
    """
    wavenumbers = np.array(spectrum['wavenumbers'])

    # Get intensities - try corrected_avg first, fall back to intensities
    if 'corrected_avg' in spectrum:
        intensities = np.array(spectrum['corrected_avg'])
    elif 'corrected' in spectrum:
        intensities = np.array(spectrum['corrected'])
    else:
        intensities = np.array(spectrum['intensities'])

    # Find indices within the normalization range
    min_wn, max_wn = wavenumber_range
    mask = (wavenumbers >= min_wn) & (wavenumbers <= max_wn)
    range_indices = np.where(mask)[0]

    if len(range_indices) == 0:
        raise ValueError(
            f"No data points found in wavenumber range {wavenumber_range}"
        )

    # Extract the region for normalization
    region = intensities[range_indices]

    # Calculate L2 norm (Euclidean norm) of the region
    l2_norm = np.linalg.norm(region)

    if l2_norm == 0:
        raise ValueError("L2 norm is zero - cannot normalize")

    # Normalize the entire spectrum by the L2 norm
    normalized = intensities / l2_norm

    return {
        'wavenumbers': wavenumbers.tolist(),
        'normalized': normalized.tolist(),
        'original': intensities.tolist(),
        'norm_factor': float(l2_norm),
        'norm_range': wavenumber_range
    }


def normalize_spectra_l2(
    sample: dict,
    references: dict,
    wavenumber_range: tuple[float, float]
) -> dict:
    """
    Normalize sample and multiple reference spectra using L2 normalization.

    Args:
        sample: Sample spectrum dictionary
        references: Dictionary of {molecule: spectrum_dict}
        wavenumber_range: Tuple of (min_wn, max_wn) for normalization

    Returns:
        Dictionary with:
        - 'sample': normalized sample spectrum
        - 'references': dict of {molecule: normalized reference spectrum}
        - 'wavenumber_range': the range used for normalization
    """
    # Normalize the sample
    normalized_sample = normalize_l2(sample, wavenumber_range)

    # Normalize each reference
    normalized_references = {}
    for molecule, ref_spectrum in references.items():
        normalized_references[molecule] = normalize_l2(ref_spectrum, wavenumber_range)

    return {
        'sample': normalized_sample,
        'references': normalized_references,
        'wavenumber_range': wavenumber_range
    }
