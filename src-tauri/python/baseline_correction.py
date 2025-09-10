#!/usr/bin/env python3
"""
baseline correction algorithm for Raman Tools.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def als_baseline(spectrum, lambda_param=1e7, p=0.01, d=2, max_iterations=10):
    """
    Asymmetric Least Squares (ALS) baseline estimation.
    Based on the algorithm described in Eilers & Boelens (2005).

    Args:
        spectrum: Input spectrum intensity values (list or array)
        lambda_param: Smoothing parameter (larger = smoother baseline)
        p: Asymmetry parameter (smaller = more weight to points below baseline)
        d: Order of differences for penalty matrix
        max_iterations: Maximum number of iterations

    Returns:
        tuple: (corrected_spectrum, baseline) as lists
    """
    spectrum = np.array(spectrum, dtype=float)
    n = len(spectrum)

    # Step 1: Construct difference matrix D (Eilers & Boelens eq. 2-3)
    # For second-order differences (d=2): D[i] = y[i] - 2*y[i+1] + y[i+2]
    if d == 2:
        D = sparse.lil_matrix((n - 2, n))
        for i in range(n - 2):
            D[i, i] = 1
            D[i, i + 1] = -2
            D[i, i + 2] = 1
        D = D.tocsr()
    else:
        # For other orders (generalized form)
        diagonals = []
        offsets = []
        for k in range(d + 1):
            from math import comb

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

    corrected_spectrum = spectrum - baseline
    return corrected_spectrum.tolist(), baseline.tolist()


def moving_average_denoise(spectrum, window_size=5):
    """
    Apply moving average denoising to spectrum.

    Args:
        spectrum: Input spectrum intensity values
        window_size: Size of the moving average window (must be odd)

    Returns:
        list: Denoised spectrum
    """
    if window_size % 2 == 0:
        window_size += 1  # Ensure odd window size

    spectrum = np.array(spectrum, dtype=float)
    half_window = window_size // 2
    denoised = np.zeros_like(spectrum)

    for i in range(len(spectrum)):
        start = max(0, i - half_window)
        end = min(len(spectrum), i + half_window + 1)
        denoised[i] = np.mean(spectrum[start:end])

    return denoised.tolist()


def process_spectrum(
    spectrum, denoise=False, window_size=5, lambda_param=1e7, p=0.01, d=2
):
    """
    Process a single spectrum with optional denoising and baseline correction.

    Args:
        spectrum: List of intensity values
        denoise: Boolean, whether to apply denoising
        window_size: Denoising window size
        lambda_param: ALS lambda parameter
        p: ALS asymmetry parameter
        d: ALS order of differences

    Returns:
        Tuple of (corrected, baseline, denoised_or_none)
    """
    # Optional denoising
    if denoise:
        denoised = moving_average_denoise(spectrum, window_size)
    else:
        denoised = spectrum

    # Baseline correction
    corrected, baseline = als_baseline(denoised, lambda_param, p, d)

    return corrected, baseline, (denoised if denoise else None)
