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
    sample_name: str,
    sample_spectrum: dict,
    reference_spectra: dict,
    molecules: list[str],
    wavenumber_range: tuple[float, float],
    output_path: str
) -> None:
    """
    Create two-panel plot showing normalization before and after.

    Top panel: Original spectra at full heights with normalization window highlighted
    Bottom panel: All normalized spectra overlaid with expected peaks marked

    Args:
        sample_name: Name of the sample (for title)
        sample_spectrum: Dictionary with 'wavenumbers', 'original', 'normalized'
        reference_spectra: Dict of {(molecule, conjugate): spectrum_dict with 'original', 'normalized'}
        molecules: List of molecules to plot (just molecule names, not tuples)
        wavenumber_range: Tuple of (min_wn, max_wn) for normalization range
        output_path: Path where the plot will be saved
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

    wavenumbers = sample_spectrum['wavenumbers']

    # ============================================================
    # Top panel: Original spectra with normalization window
    # ============================================================

    # Plot sample original (black)
    ax1.plot(wavenumbers, sample_spectrum['original'],
             color='black', linewidth=2, label=sample_name, alpha=0.8)

    # Plot reference originals (molecule colors)
    for molecule in molecules:
        # Find the reference with matching molecule (references now have tuple keys)
        for (mol, conj), ref_data in reference_spectra.items():
            if mol == molecule:
                ref_original = ref_data['original']
                color = get_color(molecule)
                ax1.plot(wavenumbers, ref_original, color=color, linewidth=1.5,
                        label=f'{molecule}-{conj} reference', alpha=0.8)
                break

    # Highlight normalization range with a rectangle
    y_min, y_max = ax1.get_ylim()
    ax1.axvspan(wavenumber_range[0], wavenumber_range[1],
               alpha=0.15, color='gray', zorder=0)

    # Add a text box indicating normalization range
    ax1.text(0.98, 0.95, f'Normalization range:\n{wavenumber_range[0]}-{wavenumber_range[1]} cm⁻¹',
            transform=ax1.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax1.set_ylabel('Intensity (Original)', fontsize=12)
    ax1.set_title('Before Normalization', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)

    # ============================================================
    # Bottom panel: Normalized spectra with peaks marked
    # ============================================================

    # Plot sample normalized (black, thicker)
    ax2.plot(wavenumbers, sample_spectrum['normalized'],
             color='black', linewidth=2.5, label=sample_name, alpha=0.8, zorder=10)

    # Plot references normalized (molecule colors)
    for molecule in molecules:
        # Find the reference with matching molecule (references now have tuple keys)
        for (mol, conj), ref_data in reference_spectra.items():
            if mol == molecule:
                ref_normalized = ref_data['normalized']
                color = get_color(molecule)

                ax2.plot(wavenumbers, ref_normalized, color=color, linewidth=1.5,
                        label=f'{molecule}-{conj} reference', alpha=0.8, zorder=5)

                # Mark expected peak
                expected_peak = get_peak(molecule)
                ax2.axvline(x=expected_peak, color=color, linestyle='--', linewidth=1.5,
                           alpha=0.5, zorder=1)
                break

    ax2.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax2.set_ylabel('Normalized Intensity', fontsize=12)
    ax2.set_title('After L2 Normalization', fontsize=12, fontweight='bold')
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(wavenumber_range[0], wavenumber_range[1])

    # Main title
    fig.suptitle(f'Normalization - {sample_name}', fontsize=14, fontweight='bold', y=0.995)

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_deconvolution(
    sample_name: str,
    sample_spectrum: dict,
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
        sample_name: Name of the sample (for title)
        sample_spectrum: Dict with 'wavenumbers' and 'normalized'
        result: Deconvolution result dict with:
                - 'contributions': {(molecule, conjugate): percentage}
                - 'reconstructed': reconstructed spectrum
                - 'residual': residual values
                - 'individual_contributions': {(molecule, conjugate): contribution}
                - 'metrics': {'rmse': value}
        wavenumber_range: Tuple of (min_wn, max_wn) for analysis range
        output_path: Path where the plot will be saved
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))

    wavenumbers = sample_spectrum['wavenumbers']
    multiplex_data = sample_spectrum['normalized']
    reconstructed = result['reconstructed']
    contributions = result['contributions']
    individual_contributions = result['individual_contributions']
    rmse = result['metrics']['rmse']

    # Get analysis range indices for residual plotting
    analysis_range = result['analysis_range']
    range_indices = analysis_range['indices']
    residual = result['residual']

    # Top panel: Multiplex vs Fitted
    ax1.plot(wavenumbers, multiplex_data, 'k-', linewidth=2,
             label='Multiplex', alpha=0.7)
    ax1.plot(wavenumbers, reconstructed, 'r--', linewidth=2,
             label='Fitted', alpha=0.8)

    # Add expected peak markers
    for mol_conj in sorted(contributions.keys()):
        molecule, conjugate = mol_conj
        expected_peak = get_peak(molecule)
        color = get_color(molecule)
        ax1.axvline(x=expected_peak, color=color, linestyle=':', linewidth=1.5,
                   alpha=0.6, zorder=1)

    ax1.set_ylabel('Intensity', fontsize=12)
    ax1.set_title(f'Multiplex Spectrum Deconvolution - {sample_name}',
                  fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(wavenumber_range[0], wavenumber_range[1])

    # Middle panel: Individual contributions
    for mol_conj in sorted(contributions.keys()):
        molecule, conjugate = mol_conj
        percentage = contributions[mol_conj]
        color = get_color(molecule)
        contribution = individual_contributions[mol_conj]

        ax2.plot(wavenumbers, contribution, color=color, linewidth=2,
                 label=f'{molecule}-{conjugate} ({percentage:.1f}%)', alpha=0.8)

        # Add expected peak marker for this molecule
        expected_peak = get_peak(molecule)
        ax2.axvline(x=expected_peak, color=color, linestyle=':', linewidth=1.5,
                   alpha=0.6, zorder=1)

    ax2.set_ylabel('Intensity', fontsize=12)
    ax2.set_title('Individual SERS Tag Contributions', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(wavenumber_range[0], wavenumber_range[1])

    # Bottom panel: Residuals (only in analysis range)
    # Create wavenumbers array for residual (only analysis range)
    residual_wavenumbers = [wavenumbers[i] for i in range_indices]

    ax3.plot(residual_wavenumbers, residual, 'g-', linewidth=1)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)

    # Add expected peak markers
    for mol_conj in sorted(contributions.keys()):
        molecule, conjugate = mol_conj
        expected_peak = get_peak(molecule)
        color = get_color(molecule)
        ax3.axvline(x=expected_peak, color=color, linestyle=':', linewidth=1.5,
                   alpha=0.6, zorder=1)

    ax3.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax3.set_ylabel('Residual', fontsize=12)
    ax3.set_title(f'Fitting Residuals (RMS: {rmse:.2f})', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(wavenumber_range[0], wavenumber_range[1])

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_deconvolution_original_scale(
    sample_name: str,
    sample_spectrum: dict,
    result: dict,
    wavenumber_range: tuple[float, float],
    output_path: str
) -> None:
    """
    Create 3-panel deconvolution plot with original-scale intensities.

    Similar to plot_deconvolution(), but shows intensities in the original scale
    (before normalization) rather than normalized scale.

    Creates a figure with three subplots:
    1. Top: Multiplex spectrum vs fitted reconstruction (original scale)
    2. Middle: Individual SERS tag contributions with percentages (original scale)
    3. Bottom: Fitting residuals with RMS value (normalized scale)

    Args:
        sample_name: Name of the sample (for title)
        sample_spectrum: Dict with 'wavenumbers', 'normalized', 'original', and 'norm_factor'
        result: Deconvolution result dict with:
                - 'contributions': {(molecule, conjugate): percentage}
                - 'reconstructed': reconstructed spectrum (normalized)
                - 'residual': residual values
                - 'individual_contributions': {(molecule, conjugate): contribution (normalized)}
                - 'metrics': {'rmse': value}
                - 'norm_factor': normalization factor
        wavenumber_range: Tuple of (min_wn, max_wn) for analysis range
        output_path: Path where the plot will be saved
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))

    wavenumbers = sample_spectrum['wavenumbers']
    original_data = sample_spectrum['original']
    contributions = result['contributions']
    norm_factor = result['norm_factor']
    rmse = result['metrics']['rmse']

    # Get analysis range indices for residual plotting
    analysis_range = result['analysis_range']
    range_indices = analysis_range['indices']
    residual = result['residual']

    # Convert normalized data back to original scale
    reconstructed_original = np.array(result['reconstructed']) * norm_factor

    # Top panel: Multiplex vs Fitted (original scale)
    ax1.plot(wavenumbers, original_data, 'k-', linewidth=2,
             label='Multiplex', alpha=0.7)
    ax1.plot(wavenumbers, reconstructed_original, 'r--', linewidth=2,
             label='Fitted', alpha=0.8)

    # Add expected peak markers
    for mol_conj in sorted(contributions.keys()):
        molecule, conjugate = mol_conj
        expected_peak = get_peak(molecule)
        color = get_color(molecule)
        ax1.axvline(x=expected_peak, color=color, linestyle=':', linewidth=1.5,
                   alpha=0.6, zorder=1)

    ax1.set_ylabel('Intensity (Original Scale)', fontsize=12)
    ax1.set_title(f'Multiplex Spectrum Deconvolution - {sample_name} (Original Scale)',
                  fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(wavenumber_range[0], wavenumber_range[1])

    # Middle panel: Individual contributions (original scale)
    individual_contributions = result['individual_contributions']

    for mol_conj in sorted(contributions.keys()):
        molecule, conjugate = mol_conj
        percentage = contributions[mol_conj]
        color = get_color(molecule)
        # Scale normalized contribution back to original
        contribution_original = np.array(individual_contributions[mol_conj]) * norm_factor

        ax2.plot(wavenumbers, contribution_original, color=color, linewidth=2,
                 label=f'{molecule}-{conjugate} ({percentage:.1f}%)', alpha=0.8)

        # Add expected peak marker for this molecule
        expected_peak = get_peak(molecule)
        ax2.axvline(x=expected_peak, color=color, linestyle=':', linewidth=1.5,
                   alpha=0.6, zorder=1)

    ax2.set_ylabel('Intensity (Original Scale)', fontsize=12)
    ax2.set_title('Individual SERS Tag Contributions', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(wavenumber_range[0], wavenumber_range[1])

    # Bottom panel: Residuals (original scale)
    # Create wavenumbers array for residual (only analysis range)
    residual_wavenumbers = [wavenumbers[i] for i in range_indices]

    # Scale residual back to original scale
    residual_original = np.array(residual) * norm_factor

    ax3.plot(residual_wavenumbers, residual_original, 'g-', linewidth=1)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)

    # Add expected peak markers
    for mol_conj in sorted(contributions.keys()):
        molecule, conjugate = mol_conj
        expected_peak = get_peak(molecule)
        color = get_color(molecule)
        ax3.axvline(x=expected_peak, color=color, linestyle=':', linewidth=1.5,
                   alpha=0.6, zorder=1)

    ax3.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12)
    ax3.set_ylabel('Residual (Original Scale)', fontsize=12)
    ax3.set_title('Fitting Residuals', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(wavenumber_range[0], wavenumber_range[1])

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_deconvolution_boxplots(
    samples: dict,
    deconv_results: dict,
    output_path: str
) -> None:
    """
    Create box plots showing distribution of molecular contributions across replicates.

    Creates a figure with box plots for each sample, showing the distribution
    of contributions for each molecule across all replicates.

    Args:
        samples: Dictionary of sample data (from load_and_process_sample)
        deconv_results: Dictionary of replicate deconvolution results
        output_path: Path where the plot will be saved
    """
    import numpy as np

    # Determine how many samples we have
    n_samples = len(samples)

    # Create subplots - one per sample
    fig, axes = plt.subplots(1, n_samples, figsize=(6 * n_samples, 6), squeeze=False)
    axes = axes.flatten()

    for idx, (sample_key, replicate_results) in enumerate(deconv_results.items()):
        ax = axes[idx]
        sample_data = samples[sample_key]
        display_name = sample_data['name']
        sample_mol_conj = sample_data['molecule_conjugates']

        # Collect contributions for each molecule-conjugate
        contributions_by_mol_conj = {mol_conj: [] for mol_conj in sample_mol_conj}

        for result in replicate_results:
            for mol_conj in sample_mol_conj:
                contributions_by_mol_conj[mol_conj].append(result['contributions'][mol_conj])

        # Prepare data for box plot
        box_data = []
        labels = []
        colors = []

        for mol_conj in sample_mol_conj:
            molecule, conjugate = mol_conj
            box_data.append(contributions_by_mol_conj[mol_conj])
            labels.append(f'{molecule}-{conjugate}')
            colors.append(get_color(molecule))

        # Create box plot
        bp = ax.boxplot(box_data, labels=labels, patch_artist=True,
                        showmeans=True, meanline=True)

        # Color the boxes
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        # Style the plot
        ax.set_ylabel('Contribution (%)', fontsize=12)
        ax.set_title(f'{display_name}\n(n={len(replicate_results)} replicates)',
                     fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_ylim(0, 100)

    # Main title
    fig.suptitle('Molecular Contributions Distribution Across Replicates',
                 fontsize=14, fontweight='bold', y=0.98)

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_peak_intensity_histogram(
    molecule: str,
    conjugate_intensities: dict[str, list[float]],
    output_path: str,
    bin_size: int = 25,
    x_max: float = None,
    y_max: int = None,
    x_label: str = "Intensity (a.u.)"
) -> None:
    """
    Create histogram of peak intensities for a molecule across conjugate types.

    Creates a histogram showing the distribution of intensities at the molecule's
    characteristic peak across all replicates, grouped by conjugate type.

    Args:
        molecule: Molecule name (e.g., "MBA", "DTNB", "TFMBA")
        conjugate_intensities: Dict mapping conjugate names to lists of intensities
                              (e.g., {"EpCAM": [120.5, 115.3, ...], "BSA": [...]})
        output_path: Full path where the plot will be saved
        bin_size: Bin size for histogram (default: 25)
        x_max: Optional maximum value for x-axis (intensity). If None, auto-scales.
        y_max: Optional maximum value for y-axis (count). If None, auto-scales.
        x_label: Label for x-axis (default: "Peak Intensity (a.u.)")
    """
    from .molecules import CONJUGATE_COLORS

    peak_wn = get_peak(molecule)

    # Create histogram plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Determine bin edges
    all_intensities = []
    for intensities in conjugate_intensities.values():
        all_intensities.extend(intensities)

    if not all_intensities:
        print(f"  Warning: No data for {molecule}, skipping")
        plt.close()
        return

    min_intensity = 0
    # Use provided x_max or calculate from data
    if x_max is not None:
        max_intensity = x_max
    else:
        max_intensity = max(all_intensities)
    bins = np.arange(min_intensity, max_intensity + bin_size, bin_size)

    # Sort conjugates alphabetically for consistent ordering
    sorted_conjugates = sorted(conjugate_intensities.keys())

    # Plot histogram for each conjugate
    for conjugate in sorted_conjugates:
        intensities = conjugate_intensities[conjugate]
        n_samples = len(intensities)

        # Create histogram
        counts, bin_edges = np.histogram(intensities, bins=bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Get color from CONJUGATE_COLORS, default to gray
        line_color = CONJUGATE_COLORS.get(conjugate, 'gray')

        # Plot as line with markers
        ax.plot(bin_centers, counts, marker='o', linewidth=2,
               label=f'{conjugate} (n={n_samples})', color=line_color)

    # Formatting
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel('Number of samples', fontsize=12)
    ax.set_title(f'{molecule} Frequency Distribution at {peak_wn} cm⁻¹',
                fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Set axis limits
    if x_max is not None:
        ax.set_xlim(0, x_max)
    else:
        ax.set_xlim(left=0)

    if y_max is not None:
        ax.set_ylim(0, y_max)
    else:
        ax.set_ylim(bottom=0)

    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
