"""
High-level workflow functions for Raman spectroscopy analysis.

These functions encapsulate common multi-step workflows to simplify experiment scripts.
"""

import os
from .io import load_spectra, ensure_output_subdir
from .baseline import apply_baseline_correction
from .averaging import calculate_average
from .plotting import plot_reference, plot_sample


def load_and_process_reference(
    directory: str,
    molecule: str,
    output_dir: str
) -> dict:
    """
    Load reference spectra, apply ALS baseline correction, average, and plot.

    This function encapsulates the complete workflow for processing reference samples:
    1. Load all spectra from the directory
    2. Apply ALS baseline correction to each spectrum
    3. Calculate the average and standard deviation
    4. Generate and save a reference plot

    Args:
        directory: Path to directory containing reference spectrum files (.txt)
        molecule: Molecule tag (e.g., "MBA", "DTNB", "TFMBA")
        output_dir: Base output directory (will create 'references/' subdirectory)

    Returns:
        dict: Averaged spectrum data with the following structure:
            {
                'wavenumbers': list[float],
                'raw_avg': list[float],
                'raw_std': list[float],
                'corrected_avg': list[float],
                'corrected_std': list[float],
                'baseline_avg': list[float],
                'count': int,
                'molecule': str
            }
    """
    # Ensure output subdirectory exists
    refs_dir = ensure_output_subdir(output_dir, "references")

    # Expand tilde in directory path
    directory = os.path.expanduser(directory)

    print(f"\nProcessing: {molecule}")
    print(f"Loading from: {directory}")

    # Load all spectra from directory
    spectra = load_spectra(directory)
    print(f"✓ Loaded {len(spectra)} spectra")

    # Apply baseline correction to each spectrum
    print("Applying baseline correction...")
    corrected_spectra = [apply_baseline_correction(s) for s in spectra]
    print("✓ Baseline correction applied")

    # Calculate average
    print("Calculating average...")
    averaged = calculate_average(corrected_spectra)
    print(f"✓ Average calculated from {averaged['count']} spectra")

    # Add molecule tag
    averaged['molecule'] = molecule

    # Generate plot
    plot_path = f"{refs_dir}/{molecule}.png"
    print(f"Creating plot: {plot_path}")
    plot_reference(averaged, molecule=molecule, output_path=plot_path)
    print("✓ Plot saved")

    return averaged


def load_and_process_sample(
    directory: str,
    name: str,
    molecules: list[str],
    output_dir: str
) -> dict:
    """
    Load sample spectra, apply ALS baseline correction, average, and plot.

    This function encapsulates the complete workflow for processing sample data:
    1. Load all spectra from the directory
    2. Apply ALS baseline correction to each spectrum
    3. Calculate the average and standard deviation
    4. Generate and save a sample plot

    Args:
        directory: Path to directory containing sample spectrum files (.txt)
        name: Display name for the sample (e.g., "Multiplex Ab 1")
        molecules: List of molecule tags present in sample (e.g., ["MBA", "DTNB", "TFMBA"])
        output_dir: Base output directory (will create 'samples/' subdirectory)

    Returns:
        dict: Averaged spectrum data with the following structure:
            {
                'wavenumbers': list[float],
                'raw_avg': list[float],
                'raw_std': list[float],
                'corrected_avg': list[float],
                'corrected_std': list[float],
                'baseline_avg': list[float],
                'count': int,
                'molecules': list[str],
                'name': str
            }
    """
    # Ensure output subdirectory exists
    samples_dir = ensure_output_subdir(output_dir, "samples")

    # Expand tilde in directory path
    directory = os.path.expanduser(directory)

    print(f"\nProcessing: {name}")
    print(f"Loading from: {directory}")

    # Load all spectra from directory
    spectra = load_spectra(directory)
    print(f"✓ Loaded {len(spectra)} spectra")

    # Apply baseline correction to each spectrum
    print("Applying baseline correction...")
    corrected_spectra = [apply_baseline_correction(s) for s in spectra]
    print("✓ Baseline correction applied")

    # Calculate average
    print("Calculating average...")
    averaged = calculate_average(corrected_spectra)
    print(f"✓ Average calculated from {averaged['count']} spectra")

    # Add metadata tags
    averaged['molecules'] = molecules
    averaged['name'] = name

    # Generate plot - create safe filename from name
    safe_name = name.replace(" ", "_").replace("/", "-")
    plot_path = f"{samples_dir}/{safe_name}.png"
    print(f"Creating plot: {plot_path}")
    plot_sample(averaged, sample_name=name, molecules=molecules, output_path=plot_path)
    print("✓ Plot saved")

    return averaged
