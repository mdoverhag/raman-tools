"""
High-level workflow functions for Raman spectroscopy analysis.

These functions encapsulate common multi-step workflows to simplify experiment scripts.
"""

import os
from .io import load_spectra, ensure_output_subdir
from .baseline import apply_baseline_correction
from .averaging import calculate_average
from .plotting import plot_reference, plot_sample, plot_normalization, plot_deconvolution, plot_deconvolution_boxplots
from .normalization import normalize_spectra_l2
from .deconvolution import deconvolve_nnls


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
                'name': str,
                'replicates': list[dict]  # List of individual corrected spectra
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

    # Add metadata tags to averaged data
    averaged['molecules'] = molecules
    averaged['name'] = name

    # Add replicates to the result
    averaged['replicates'] = corrected_spectra

    # Generate plot - create safe filename from name
    safe_name = name.replace(" ", "_").replace("/", "-")
    plot_path = f"{samples_dir}/{safe_name}.png"
    print(f"Creating plot: {plot_path}")
    plot_sample(averaged, sample_name=name, molecules=molecules, output_path=plot_path)
    print("✓ Plot saved")

    return averaged


def normalize_and_deconvolve_samples(
    samples: dict,
    references: dict,
    wavenumber_range: tuple[float, float],
    output_dir: str
) -> dict:
    """
    Normalize and deconvolve all samples against references.

    This function encapsulates the complete workflow for normalization and deconvolution:
    1. Normalize each sample and all references using L2 normalization
    2. Create normalization plot for each sample
    3. Perform NNLS deconvolution on each normalized sample
    4. Create deconvolution plot for each sample
    5. Print progress and results for each sample

    Args:
        samples: Dictionary of sample data (from load_and_process_sample)
        references: Dictionary of reference data (from load_and_process_reference)
        wavenumber_range: Tuple of (min, max) wavenumber for analysis window
        output_dir: Base output directory (will create subdirectories)

    Returns:
        dict: Deconvolution results for each sample (list of results for all replicates):
            {
                sample_key: [
                    {
                        'contributions': {molecule: percentage},
                        'reconstructed': [...],
                        'residual': [...],
                        'metrics': {'rmse': float, 'r_squared': float}
                    },
                    # ... one result per replicate
                ]
            }
    """
    # Validate all required references exist before processing
    for sample_key, sample_data in samples.items():
        for molecule in sample_data['molecules']:
            if molecule not in references:
                raise ValueError(
                    f"Sample '{sample_data['name']}' requires molecule '{molecule}' "
                    f"but no reference was loaded for it. "
                    f"Available references: {list(references.keys())}"
                )

    # Ensure output subdirectories exist
    norm_dir = ensure_output_subdir(output_dir, "normalization")
    deconv_dir = ensure_output_subdir(output_dir, "deconvolution")

    # Store results for summary
    deconv_results = {}

    for sample_key, sample_data in samples.items():
        display_name = sample_data['name']
        sample_molecules = sample_data['molecules']
        replicates = sample_data['replicates']

        print(f"\n{display_name}:")
        print(f"  Processing {len(replicates)} replicates...")
        print(f"  Normalizing and deconvolving (range: {wavenumber_range[0]}-{wavenumber_range[1]} cm⁻¹)...")

        # Process each replicate
        replicate_results = []

        for i, replicate in enumerate(replicates):
            # Normalize this replicate
            normalized = normalize_spectra_l2(
                sample=replicate,
                references=references,
                wavenumber_range=wavenumber_range
            )

            # Deconvolve this replicate
            deconv_result = deconvolve_nnls(
                sample_spectrum=normalized['sample'],
                reference_spectra=normalized['references'],
                wavenumber_range=wavenumber_range
            )

            replicate_results.append(deconv_result)

        # Store all replicate results
        deconv_results[sample_key] = replicate_results

        # Create plots using the averaged spectrum (for backward compatibility)
        # Normalize the averaged spectrum for plotting
        normalized_avg = normalize_spectra_l2(
            sample=sample_data,
            references=references,
            wavenumber_range=wavenumber_range
        )

        plot_normalization(
            sample_name=f"{display_name} (averaged)",
            sample_spectrum=normalized_avg['sample'],
            reference_spectra=normalized_avg['references'],
            molecules=sample_molecules,
            wavenumber_range=wavenumber_range,
            output_path=f"{norm_dir}/{sample_key}_averaged.png"
        )

        # Deconvolve the averaged spectrum for plotting
        deconv_avg = deconvolve_nnls(
            sample_spectrum=normalized_avg['sample'],
            reference_spectra=normalized_avg['references'],
            wavenumber_range=wavenumber_range
        )

        plot_deconvolution(
            sample_name=f"{display_name} (averaged)",
            sample_spectrum=normalized_avg['sample'],
            result=deconv_avg,
            wavenumber_range=wavenumber_range,
            output_path=f"{deconv_dir}/{sample_key}_averaged.png"
        )

        # Calculate and print summary statistics across replicates
        import numpy as np

        contributions_by_mol = {mol: [] for mol in sample_molecules}
        r_squared_values = []

        for result in replicate_results:
            for mol in sample_molecules:
                contributions_by_mol[mol].append(result['contributions'][mol])
            r_squared_values.append(result['metrics']['r_squared'])

        # Print mean ± std for each molecule
        contributions_str = ", ".join([
            f"{mol}={np.mean(contributions_by_mol[mol]):.1f}±{np.std(contributions_by_mol[mol]):.1f}%"
            for mol in sample_molecules
        ])
        r2_mean = np.mean(r_squared_values)
        r2_std = np.std(r_squared_values)

        print(f"  ✓ Processed {len(replicates)} replicates")
        print(f"  ✓ Contributions: {contributions_str}")
        print(f"  ✓ R²: {r2_mean:.3f}±{r2_std:.3f}")

    # Create box plots showing distribution across all samples
    print(f"\nCreating box plots...")
    boxplot_path = f"{deconv_dir}/boxplots.png"
    plot_deconvolution_boxplots(
        samples=samples,
        deconv_results=deconv_results,
        output_path=boxplot_path
    )
    print(f"✓ Box plots saved: {boxplot_path}")

    return deconv_results


def print_experiment_summary(
    output_dir: str,
    references: dict,
    samples: dict,
    deconv_results: dict
) -> None:
    """
    Print a summary of the experiment results.

    This function prints a formatted summary including:
    1. Experiment completion header
    2. Output directory location
    3. References processed (molecule names and spectrum counts)
    4. Samples processed (sample names and spectrum counts)
    5. Deconvolution results table (contributions and R² values)
    6. Plots saved locations

    Args:
        output_dir: Path to the output directory
        references: Dictionary of reference data (from load_and_process_reference)
        samples: Dictionary of sample data (from load_and_process_sample)
        deconv_results: Dictionary of deconvolution results (from normalize_and_deconvolve_samples)
    """
    print(f"\nOutput directory: {output_dir}")

    print(f"\nReferences processed:")
    for molecule, data in references.items():
        print(f"  {molecule}: {data['count']} spectra")

    print(f"\nSamples processed:")
    for sample_key, data in samples.items():
        print(f"  {data['name']}: {data['count']} spectra")

    # Determine all unique molecules across all samples
    all_molecules = set()
    for sample_data in samples.values():
        all_molecules.update(sample_data['molecules'])
    all_molecules = sorted(all_molecules)

    print(f"\nDeconvolution Results (mean ± std across replicates):")

    # Build header dynamically - need more space for "mean±std" format
    header = f"{'Sample':<25}"
    for mol in all_molecules:
        header += f" {mol:>16}"
    header += f" {'R²':>14}"
    print(f"\n{header}")
    print("-" * (25 + len(all_molecules) * 17 + 15))

    # Print each sample's results (now with mean ± std)
    import numpy as np

    for sample_key, replicate_results in deconv_results.items():
        display_name = samples[sample_key]['name']
        sample_molecules = samples[sample_key]['molecules']

        # Calculate mean ± std for each molecule
        contributions_by_mol = {mol: [] for mol in sample_molecules}
        r_squared_values = []

        for result in replicate_results:
            for mol in sample_molecules:
                contributions_by_mol[mol].append(result['contributions'][mol])
            r_squared_values.append(result['metrics']['r_squared'])

        # Build row with mean±std - fixed width format
        row = f"{display_name:<25}"
        for mol in all_molecules:
            if mol in sample_molecules:
                mean = np.mean(contributions_by_mol[mol])
                std = np.std(contributions_by_mol[mol])
                # Fixed width: "XX.X±XX.X%" = 11 chars, right-aligned in 16 char field
                row += f" {f'{mean:.1f}±{std:.1f}%':>16}"
            else:
                row += f" {'-':>16}"

        r2_mean = np.mean(r_squared_values)
        r2_std = np.std(r_squared_values)
        # Fixed width for R²
        row += f" {f'{r2_mean:.3f}±{r2_std:.3f}':>14}"
        print(row)

    print(f"\nPlots saved in:")
    print(f"  {output_dir}/references/")
    print(f"  {output_dir}/samples/")
    print(f"  {output_dir}/normalization/")
    print(f"  {output_dir}/deconvolution/")
