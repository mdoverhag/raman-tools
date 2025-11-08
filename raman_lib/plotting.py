"""
Plotting functions for Raman spectroscopy data.

Creates publication-quality plots for references, samples, normalization,
and deconvolution results.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from .molecules import get_peak, get_color


def plot_reference(
    averaged_spectrum: dict,
    molecule: str,
    output_path: str
) -> None:
    """
    Plot a reference spectrum showing raw, corrected, and expected peak.

    Creates a plot with:
    - Average raw spectrum (light gray)
    - Average corrected spectrum (molecule color)
    - Vertical line at expected peak for the molecule

    Args:
        averaged_spectrum: Dictionary with 'wavenumbers', 'raw_avg', 'corrected_avg'
        molecule: Molecule name (e.g., "MBA", "DTNB", "TFMBA")
        output_path: Path where the plot will be saved
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    wavenumbers = averaged_spectrum['wavenumbers']
    raw_avg = averaged_spectrum['raw_avg']
    corrected_avg = averaged_spectrum['corrected_avg']

    # Get molecule properties
    expected_peak = get_peak(molecule)
    color = get_color(molecule)

    # Plot raw spectrum (lighter, in background)
    ax.plot(wavenumbers, raw_avg, color='lightgray', linewidth=1.5,
            label='Raw (average)', alpha=0.8)

    # Plot corrected spectrum (molecule color, prominent)
    ax.plot(wavenumbers, corrected_avg, color=color, linewidth=2,
            label='Baseline corrected (average)')

    # Mark expected peak with vertical line
    ax.axvline(x=expected_peak, color=color, linestyle='--', linewidth=1.5,
               alpha=0.7, label=f'Expected peak ({expected_peak} cm⁻¹)')

    # Labels and formatting
    ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Intensity', fontsize=12)
    ax.set_title(f'Reference Spectrum - {molecule}', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_sample(
    averaged_spectrum: dict,
    sample_name: str,
    molecules: list[str],
    output_path: str
) -> None:
    """
    Plot a sample spectrum showing raw, corrected, and expected peaks.

    Creates a plot with:
    - Average raw spectrum (light gray)
    - Average corrected spectrum (black)
    - Vertical lines at expected peaks for all molecules in the sample

    Args:
        averaged_spectrum: Dictionary with 'wavenumbers', 'raw_avg', 'corrected_avg'
        sample_name: Name of the sample (e.g., "Multi Ab 1")
        molecules: List of molecules present (e.g., ["MBA", "DTNB", "TFMBA"])
        output_path: Path where the plot will be saved
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    wavenumbers = averaged_spectrum['wavenumbers']
    raw_avg = averaged_spectrum['raw_avg']
    corrected_avg = averaged_spectrum['corrected_avg']

    # Plot raw spectrum (lighter, in background)
    ax.plot(wavenumbers, raw_avg, color='lightgray', linewidth=1.5,
            label='Raw (average)', alpha=0.8)

    # Plot corrected spectrum (black, prominent)
    ax.plot(wavenumbers, corrected_avg, color='black', linewidth=2,
            label='Baseline corrected (average)')

    # Mark expected peaks for all molecules
    for molecule in molecules:
        expected_peak = get_peak(molecule)
        color = get_color(molecule)

        ax.axvline(x=expected_peak, color=color, linestyle='--', linewidth=1.5,
                   alpha=0.7, label=f'{molecule} ({expected_peak} cm⁻¹)')

    # Labels and formatting
    ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Intensity', fontsize=12)
    ax.set_title(f'Sample Spectrum - {sample_name}', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_normalization(
    sample_spectrum: dict,
    reference_spectra: dict,
    molecules: list[str],
    wavenumber_range: tuple[float, float],
    output_path: str
) -> None:
    """
    Plot normalized sample and reference spectra.

    Creates a plot showing:
    - Normalized sample spectrum (black, thick)
    - Normalized reference spectra (molecule colors, thin)
    - Shaded region for normalization range
    - Expected peaks marked

    Args:
        sample_spectrum: Dictionary with 'wavenumbers' and 'normalized'
        reference_spectra: Dict of {molecule: spectrum_dict with 'normalized'}
        molecules: List of molecules to plot
        wavenumber_range: Tuple of (min_wn, max_wn) for normalization range
        output_path: Path where the plot will be saved
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    wavenumbers = sample_spectrum['wavenumbers']
    sample_normalized = sample_spectrum['normalized']

    # Plot sample (black, prominent)
    ax.plot(wavenumbers, sample_normalized, color='black', linewidth=2.5,
            label='Sample (normalized)', zorder=10)

    # Plot references (molecule colors)
    for molecule in molecules:
        if molecule in reference_spectra:
            ref_normalized = reference_spectra[molecule]['normalized']
            color = get_color(molecule)

            ax.plot(wavenumbers, ref_normalized, color=color, linewidth=1.5,
                    label=f'{molecule} reference', alpha=0.8, zorder=5)

            # Mark expected peak
            expected_peak = get_peak(molecule)
            ax.axvline(x=expected_peak, color=color, linestyle=':', linewidth=1,
                       alpha=0.5, zorder=1)

    # Highlight normalization range
    ax.axvspan(wavenumber_range[0], wavenumber_range[1],
               alpha=0.1, color='gray', zorder=0,
               label=f'Normalization range ({wavenumber_range[0]}-{wavenumber_range[1]} cm⁻¹)')

    # Labels and formatting
    ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax.set_ylabel('Normalized Intensity', fontsize=12)
    ax.set_title('Normalized Spectra', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_deconvolution(
    multiplex_spectrum: dict,
    reference_spectra: dict,
    result: dict,
    wavenumber_range: tuple[float, float],
    output_path: str
) -> None:
    """
    Create 3-panel deconvolution plot.

    Creates a figure with three subplots:
    1. Top: Multiplex spectrum vs fitted reconstruction
    2. Middle: Individual SERS tag contributions with percentages
    3. Bottom: Fitting residuals with RMS value

    Args:
        multiplex_spectrum: Dict with 'wavenumbers' and 'normalized'
        reference_spectra: Dict of {molecule: spectrum_dict}
        result: Deconvolution result dict with:
                - 'contributions': {molecule: percentage}
                - 'reconstructed': reconstructed spectrum
                - 'residual': residual values
                - 'metrics': {'rmse': value, 'r_squared': value}
        wavenumber_range: Tuple of (min_wn, max_wn) for analysis range
        output_path: Path where the plot will be saved
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))

    wavenumbers = multiplex_spectrum['wavenumbers']
    multiplex_data = multiplex_spectrum['normalized']
    reconstructed = result['reconstructed']
    residual = result['residual']
    contributions = result['contributions']
    rmse = result['metrics']['rmse']

    # Top panel: Multiplex vs Fitted
    ax1.plot(wavenumbers, multiplex_data, 'k-', linewidth=2,
             label='Multiplex', alpha=0.7)
    ax1.plot(wavenumbers, reconstructed, 'r--', linewidth=2,
             label='Fitted', alpha=0.8)
    ax1.set_ylabel('Intensity', fontsize=12)
    ax1.set_title('Multiplex Spectrum Deconvolution', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Middle panel: Individual contributions
    for molecule, percentage in contributions.items():
        color = get_color(molecule)
        # Get the contribution for this molecule
        # This is the reference spectrum scaled by its coefficient
        contribution = result.get(f'{molecule}_contribution', [0] * len(wavenumbers))

        ax2.plot(wavenumbers, contribution, color=color, linewidth=2,
                 label=f'{molecule} ({percentage:.1f}%)', alpha=0.8)

    ax2.set_ylabel('Intensity', fontsize=12)
    ax2.set_title('Individual SERS Tag Contributions', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Bottom panel: Residuals
    ax3.plot(wavenumbers, residual, 'g-', linewidth=1)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax3.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax3.set_ylabel('Residual', fontsize=12)
    ax3.set_title(f'Fitting Residuals (RMS: {rmse:.2f})', fontsize=12)
    ax3.grid(True, alpha=0.3)

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
