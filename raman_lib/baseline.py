"""
Baseline correction algorithms for Raman spectroscopy.

Implements the ALS (Asymmetric Least Squares) baseline correction method
based on Eilers & Boelens (2005).
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def als_baseline(
    intensities: list[float],
    lambda_param: float = 1e7,
    p: float = 0.01,
    d: int = 2,
    max_iterations: int = 10
) -> tuple[list[float], list[float]]:
    """
    Apply Asymmetric Least Squares (ALS) baseline correction.

    This algorithm iteratively fits a smooth baseline to the spectrum by
    penalizing deviations from smoothness and asymmetrically weighting points
    above (peaks) and below (baseline) the fitted curve.

    Based on: Eilers & Boelens (2005), "Baseline Correction with Asymmetric
    Least Squares Smoothing"

    Args:
        intensities: List of intensity values
        lambda_param: Smoothing parameter (larger = smoother baseline).
                     Typical range: 1e5 to 1e9. Default: 1e7
        p: Asymmetry parameter (smaller = more weight to points below baseline).
           Typical range: 0.001 to 0.1. Default: 0.01
        d: Order of differences for penalty matrix. Default: 2 (second-order)
        max_iterations: Maximum number of iterations. Default: 10

    Returns:
        Tuple of (corrected_intensities, baseline) as lists

    Example:
        >>> corrected, baseline = als_baseline(spectrum['intensities'])
        >>> # Now spectrum peaks are isolated from the baseline
    """
    spectrum = np.array(intensities, dtype=float)
    n = len(spectrum)

    # Step 1: Construct difference matrix D (Eilers & Boelens eq. 2-3)
    # For second-order differences (d=2): D[i] = y[i] - 2*y[i+1] + y[i+2]
    if d == 2:
        # Build sparse matrix efficiently for second-order differences
        D = sparse.lil_matrix((n - 2, n))
        for i in range(n - 2):
            D[i, i] = 1
            D[i, i + 1] = -2
            D[i, i + 2] = 1
        D = D.tocsr()
    else:
        # Generalized form for other orders
        from math import comb
        diagonals = []
        offsets = []
        for k in range(d + 1):
            coeff = (-1) ** (d - k) * comb(d, k)
            diagonals.append(np.full(n - k, coeff))
            offsets.append(k)
        D = sparse.diags(diagonals, offsets, shape=(n - d, n), format="csr")

    # Step 2: Construct penalty matrix P (Eilers & Boelens eq. 4)
    # P = λ * D'D where λ controls the smoothness
    P = lambda_param * (D.T @ D)

    # Step 3: Initialize weights vector w (all ones initially)
    w = np.ones(n)

    # Step 4: Iterative reweighting (Eilers & Boelens Section 4)
    for iteration in range(max_iterations):
        # Construct diagonal weight matrix W
        W = sparse.diags(w, format="csr")

        # Solve weighted penalized least squares (Eilers & Boelens eq. 5)
        # (W + P) * z = W * y, where z is the baseline and y is the spectrum
        A = W + P
        baseline = spsolve(A, W @ spectrum)

        # Update weights asymmetrically (Eilers & Boelens eq. 6-7)
        # Points above baseline (peaks) get small weight p
        # Points below baseline get large weight (1-p)
        residuals = spectrum - baseline
        w = p * (residuals > 0) + (1 - p) * (residuals <= 0)

    # Calculate corrected spectrum by subtracting baseline
    corrected_spectrum = spectrum - baseline

    return corrected_spectrum.tolist(), baseline.tolist()


def apply_baseline_correction(
    spectrum: dict,
    lambda_param: float = 1e7,
    p: float = 0.01,
    d: int = 2
) -> dict:
    """
    Apply ALS baseline correction to a complete spectrum.

    Args:
        spectrum: Dictionary with 'wavenumbers' and 'intensities'
        lambda_param: ALS smoothing parameter
        p: ALS asymmetry parameter
        d: ALS order of differences

    Returns:
        Dictionary with corrected spectrum and baseline:
        {
            'wavenumbers': [...],
            'intensities': [...],  # original
            'corrected': [...],     # baseline-corrected
            'baseline': [...]       # extracted baseline
        }
    """
    corrected, baseline = als_baseline(
        spectrum['intensities'],
        lambda_param=lambda_param,
        p=p,
        d=d
    )

    return {
        'wavenumbers': spectrum['wavenumbers'],
        'intensities': spectrum['intensities'],
        'corrected': corrected,
        'baseline': baseline
    }
