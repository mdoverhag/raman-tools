#!/usr/bin/env python3
"""
Perform NNLS (Non-Negative Least Squares) deconvolution on multiplex SERS spectra.
"""

import json
import sys
import numpy as np
from scipy.optimize import nnls


def perform_nnls_deconvolution(multiplex_spectrum, reference_spectra,
                               wavenumber_range=(1000, 1500), total_points=1801):
    """
    Perform NNLS deconvolution to find the contribution of each reference spectrum
    to the multiplex spectrum.

    Args:
        multiplex_spectrum: Normalized multiplex spectrum (full length)
        reference_spectra: Dict of normalized reference spectra
        wavenumber_range: Tuple of (start_wavenumber, end_wavenumber)
        total_points: Total number of data points

    Returns:
        Dict containing contributions and reconstruction metrics
    """
    # Find the indices for the wavenumber range
    start_wn, end_wn = wavenumber_range
    default_start, default_end = 200.0, 2000.0

    # Calculate step size and indices
    step = (default_end - default_start) / (total_points - 1)
    start_idx = int((start_wn - default_start) / step)
    end_idx = int((end_wn - default_start) / step)

    # Clamp to valid range
    start_idx = max(0, min(start_idx, total_points - 1))
    end_idx = max(0, min(end_idx, total_points - 1))

    # Extract the region for NNLS
    multiplex_region = np.array(multiplex_spectrum[start_idx:end_idx + 1])

    # Build the reference matrix (each column is a reference spectrum)
    reference_names = list(reference_spectra.keys())
    reference_matrix = np.column_stack([
        np.array(reference_spectra[name][start_idx:end_idx + 1])
        for name in reference_names
    ])

    # Perform NNLS: find x such that ||Ax - b||^2 is minimized, x >= 0
    # where A is reference_matrix, b is multiplex_region
    coefficients, residual_norm = nnls(reference_matrix, multiplex_region)

    # Normalize coefficients to sum to 1 (get percentages)
    total = np.sum(coefficients)
    if total > 0:
        percentages = coefficients / total * 100
    else:
        percentages = coefficients

    # Reconstruct the spectrum using the found coefficients
    reconstructed = np.dot(reference_matrix, coefficients)

    # Calculate residual (difference between original and reconstructed)
    residual = multiplex_region - reconstructed

    # Calculate RMSE (Root Mean Square Error)
    rmse = np.sqrt(np.mean(residual**2))

    # Calculate R-squared (coefficient of determination)
    ss_res = np.sum(residual**2)
    ss_tot = np.sum((multiplex_region - np.mean(multiplex_region))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # Create full-length reconstructed spectrum (for visualization)
    full_reconstructed = np.zeros(total_points)
    full_reconstructed[start_idx:end_idx + 1] = reconstructed

    # Create result dictionary
    result = {
        "contributions": {
            name: float(percentages[i])
            for i, name in enumerate(reference_names)
        },
        "coefficients": {
            name: float(coefficients[i])
            for i, name in enumerate(reference_names)
        },
        "reconstructedSpectrum": full_reconstructed.tolist(),
        "residual": residual.tolist(),
        "metrics": {
            "rmse": float(rmse),
            "rSquared": float(r_squared),
            "residualNorm": float(residual_norm),
            "totalContribution": float(np.sum(percentages))
        },
        "analysisRange": {
            "startIndex": start_idx,
            "endIndex": end_idx,
            "wavenumberRange": wavenumber_range
        }
    }

    return result


def main():
    # Read input from stdin
    input_data = json.loads(sys.stdin.read())

    multiplex_spectrum = input_data["multiplexSpectrum"]
    reference_spectra = input_data["referenceSpectra"]
    wavenumber_range = tuple(input_data.get("wavenumberRange", [1000, 1500]))
    total_points = input_data.get("totalPoints", 1801)

    try:
        # Perform deconvolution
        result = perform_nnls_deconvolution(
            multiplex_spectrum,
            reference_spectra,
            wavenumber_range,
            total_points
        )

        # Output result
        print(json.dumps(result))
    except Exception as e:
        # Output error
        error_result = {
            "error": str(e)
        }
        print(json.dumps(error_result))
        sys.exit(1)


if __name__ == "__main__":
    main()