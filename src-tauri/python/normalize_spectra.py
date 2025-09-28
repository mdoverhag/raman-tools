#!/usr/bin/env python3
"""
Normalize spectra using L2 normalization within a specified wavenumber range.
"""

import json
import sys
import numpy as np


def l2_normalize(spectrum, start_idx, end_idx):
    """
    Apply L2 normalization to a spectrum within the specified index range.

    Args:
        spectrum: Array of intensity values
        start_idx: Starting index for normalization range
        end_idx: Ending index for normalization range (inclusive)

    Returns:
        Normalized spectrum (full length, but only normalized within range)
    """
    spectrum_array = np.array(spectrum)

    # Extract the region for normalization
    region = spectrum_array[start_idx : end_idx + 1]

    # Calculate L2 norm of the region
    l2_norm = np.linalg.norm(region)

    # Avoid division by zero
    if l2_norm == 0:
        return spectrum_array.tolist()

    # Normalize only the specified region
    normalized = spectrum_array.copy()
    normalized[start_idx : end_idx + 1] = region / l2_norm

    return normalized.tolist()


def find_wavenumber_indices(
    wavenumber_range, total_points=1801, default_start=200.0, default_end=2000.0
):
    """
    Find the indices corresponding to the wavenumber range.

    Args:
        wavenumber_range: Tuple of (start_wavenumber, end_wavenumber)
        total_points: Total number of data points (default 1801)
        default_start: Default starting wavenumber (200 cm⁻¹)
        default_end: Default ending wavenumber (2000 cm⁻¹)

    Returns:
        Tuple of (start_index, end_index)
    """
    start_wn, end_wn = wavenumber_range

    # Calculate step size
    step = (default_end - default_start) / (total_points - 1)

    # Calculate indices
    start_idx = int((start_wn - default_start) / step)
    end_idx = int((end_wn - default_start) / step)

    # Clamp to valid range
    start_idx = max(0, min(start_idx, total_points - 1))
    end_idx = max(0, min(end_idx, total_points - 1))

    return start_idx, end_idx


def main():
    # Read input from stdin
    input_data = json.loads(sys.stdin.read())

    multiplex_spectrum = input_data["multiplexSpectrum"]
    reference_spectra = input_data["referenceSpectra"]
    wavenumber_range = tuple(input_data["wavenumberRange"])

    # Find the indices for the wavenumber range
    start_idx, end_idx = find_wavenumber_indices(wavenumber_range)

    # Normalize the multiplex spectrum
    normalized_multiplex = l2_normalize(multiplex_spectrum, start_idx, end_idx)

    # Normalize each reference spectrum
    normalized_references = {}
    for name, spectrum in reference_spectra.items():
        normalized_references[name] = l2_normalize(spectrum, start_idx, end_idx)

    # Return results
    result = {
        "normalizedMultiplex": normalized_multiplex,
        "normalizedReferences": normalized_references,
        "normalizationRange": {
            "startIndex": start_idx,
            "endIndex": end_idx,
            "wavenumberRange": wavenumber_range,
        },
    }

    print(json.dumps(result))


if __name__ == "__main__":
    main()
