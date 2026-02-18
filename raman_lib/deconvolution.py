"""
Deconvolution of multiplex SERS spectra using NNLS.

Uses Non-Negative Least Squares to determine the contribution of each
reference spectrum to a multiplex spectrum.
"""

import numpy as np
from scipy.optimize import nnls

from .molecules import get_peak


def deconvolve_nnls(
    sample_spectrum: dict,
    reference_spectra: dict,
    wavenumber_range: tuple[float, float],
) -> dict:
    """
    Perform NNLS deconvolution on a normalized multiplex spectrum.

    Args:
        sample_spectrum: Normalized sample spectrum dict with 'wavenumbers', 'normalized', 'norm_factor'
        reference_spectra: Dict of {(molecule, conjugate): normalized_spectrum}
        wavenumber_range: Tuple of (min_wn, max_wn) for analysis range

    Returns:
        Dictionary with:
        - 'contributions': {(molecule, conjugate): percentage}
        - 'coefficients': {(molecule, conjugate): raw coefficient}
        - 'reconstructed': full-length reconstructed spectrum
        - 'residual': residual in analysis range
        - 'metrics': {'rmse': value, 'r_squared': value}
        - 'individual_contributions': {(molecule, conjugate): full-length contribution}
        - 'norm_factor': norm factor used for sample normalization
    """
    wavenumbers = np.array(sample_spectrum["wavenumbers"])
    sample_normalized = np.array(sample_spectrum["normalized"])

    # Find indices within the analysis range
    min_wn, max_wn = wavenumber_range
    mask = (wavenumbers >= min_wn) & (wavenumbers <= max_wn)
    range_indices = np.where(mask)[0]

    if len(range_indices) == 0:
        raise ValueError(f"No data points found in wavenumber range {wavenumber_range}")

    # Extract the region for NNLS
    sample_region = sample_normalized[range_indices]

    # Build the reference matrix (each column is a reference spectrum)
    # reference_keys is a list of (molecule, conjugate) tuples
    reference_keys = list(reference_spectra.keys())
    reference_matrix = np.column_stack(
        [
            np.array(reference_spectra[key]["normalized"])[range_indices]
            for key in reference_keys
        ]
    )

    # Perform NNLS: find x such that ||Ax - b||^2 is minimized, x >= 0
    # where A is reference_matrix, b is sample_region
    coefficients, residual_norm = nnls(reference_matrix, sample_region)

    # Normalize coefficients to sum to 1 (get percentages)
    total = np.sum(coefficients)
    if total > 0:
        percentages = coefficients / total * 100
    else:
        percentages = coefficients

    # Reconstruct the spectrum using the found coefficients
    reconstructed_region = np.dot(reference_matrix, coefficients)

    # Calculate residual (difference between original and reconstructed)
    residual = sample_region - reconstructed_region

    # Calculate RMSE (Root Mean Square Error)
    rmse = np.sqrt(np.mean(residual**2))

    # Calculate R-squared (coefficient of determination)
    ss_res = np.sum(residual**2)
    ss_tot = np.sum((sample_region - np.mean(sample_region)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # Create full-length reconstructed spectrum (zero outside analysis range)
    full_reconstructed = np.zeros(len(wavenumbers))
    full_reconstructed[range_indices] = reconstructed_region

    # Create individual contributions for plotting
    individual_contributions = {}
    for i, key in enumerate(reference_keys):
        # Each molecule-conjugate's contribution is its reference spectrum scaled by its coefficient
        full_contribution = np.zeros(len(wavenumbers))
        ref_full = np.array(reference_spectra[key]["normalized"])
        # Scale the entire reference spectrum by its coefficient
        full_contribution = ref_full * coefficients[i]
        individual_contributions[key] = full_contribution.tolist()

    # Get norm_factor from sample spectrum (if available)
    norm_factor = sample_spectrum.get("norm_factor", None)

    # Create result dictionary
    result = {
        "contributions": {
            key: float(percentages[i]) for i, key in enumerate(reference_keys)
        },
        "coefficients": {
            key: float(coefficients[i]) for i, key in enumerate(reference_keys)
        },
        "reconstructed": full_reconstructed.tolist(),
        "residual": residual.tolist(),
        "metrics": {
            "rmse": float(rmse),
            "r_squared": float(r_squared),
            "residual_norm": float(residual_norm),
        },
        "individual_contributions": individual_contributions,
        "analysis_range": {
            "wavenumber_range": wavenumber_range,
            "indices": range_indices.tolist(),
        },
        "norm_factor": norm_factor,
        "wavenumbers": sample_spectrum["wavenumbers"],
    }

    return result


def scale_contributions(deconv_result: dict, references: dict) -> dict:
    """
    Scale deconvolution contributions by reference signal strength.

    Different Raman reporters have very different inherent signal strengths.
    This correction divides each molecule's deconvolved signal by its reference
    peak intensity, then re-normalizes to get molecular abundance percentages.

    Args:
        deconv_result: Single replicate deconvolution result from deconvolve_nnls()
        references: Dict of {(molecule, conjugate): averaged_reference} with
                    'wavenumbers' and 'corrected_avg' keys (original, non-normalized)

    Returns:
        Dictionary with:
        - 'scaled_contributions': {(molecule, conjugate): percentage}
        - 'reference_peaks': {(molecule, conjugate): peak_intensity}
        - 'multiplex_peaks': {(molecule, conjugate): peak_intensity}
        - 'relative_amounts': {(molecule, conjugate): relative_amount}
    """
    wavenumbers = np.array(deconv_result["wavenumbers"])
    norm_factor = deconv_result["norm_factor"]

    reference_peaks = {}
    multiplex_peaks = {}
    relative_amounts = {}

    for mol_conj in deconv_result["contributions"]:
        molecule, conjugate = mol_conj
        peak_wn = get_peak(molecule)
        peak_idx = np.argmin(np.abs(wavenumbers - peak_wn))

        # Reference peak: intensity at characteristic wavenumber in averaged corrected reference
        ref = references[mol_conj]
        ref_wavenumbers = np.array(ref["wavenumbers"])
        ref_peak_idx = np.argmin(np.abs(ref_wavenumbers - peak_wn))
        reference_peaks[mol_conj] = float(np.array(ref["corrected_avg"])[ref_peak_idx])

        # Multiplex peak: intensity at characteristic wavenumber in de-normalized individual contribution
        contribution_normalized = np.array(deconv_result["individual_contributions"][mol_conj])
        contribution_original = contribution_normalized * norm_factor
        multiplex_peaks[mol_conj] = float(contribution_original[peak_idx])

        # Relative amount: multiplex peak / reference peak
        if reference_peaks[mol_conj] > 0:
            relative_amounts[mol_conj] = multiplex_peaks[mol_conj] / reference_peaks[mol_conj]
        else:
            relative_amounts[mol_conj] = 0.0

    # Normalize relative amounts to percentages
    total = sum(relative_amounts.values())
    if total > 0:
        scaled_contributions = {
            mc: ra / total * 100 for mc, ra in relative_amounts.items()
        }
    else:
        scaled_contributions = {mc: 0.0 for mc in relative_amounts}

    return {
        "scaled_contributions": scaled_contributions,
        "reference_peaks": reference_peaks,
        "multiplex_peaks": multiplex_peaks,
        "relative_amounts": relative_amounts,
    }
