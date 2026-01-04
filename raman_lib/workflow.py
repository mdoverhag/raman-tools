"""
High-level workflow functions for Raman spectroscopy analysis.

These functions encapsulate common multi-step workflows to simplify experiment scripts.
"""

import os
from .io import load_spectra, ensure_output_subdir
from .baseline import apply_baseline_correction
from .averaging import calculate_average
from .plotting import plot_reference, plot_sample, plot_normalization, plot_deconvolution, plot_deconvolution_original_scale, plot_deconvolution_boxplots
from .normalization import normalize_spectra_l2
from .deconvolution import deconvolve_nnls


def build_reference_dict(reference_list: list[dict]) -> dict:
    """
    Build a references dict from a list of reference data.

    Uses the (molecule, conjugate) tuple from each reference as the dict key.
    Validates that all molecule-conjugate pairs are unique within the reference set.

    Args:
        reference_list: List of reference dicts returned from load_and_process_reference()

    Returns:
        dict: Dictionary with (molecule, conjugate) tuples as keys and reference data as values

    Raises:
        ValueError: If duplicate molecule-conjugate pairs are found in the reference list

    Example:
        >>> references = build_reference_dict([
        ...     load_and_process_reference(..., molecule="MBA", conjugate="EpCAM", ...),
        ...     load_and_process_reference(..., molecule="DTNB", conjugate="HER2", ...),
        ... ])
        >>> # references = {("MBA", "EpCAM"): {...}, ("DTNB", "HER2"): {...}}
    """
    references = {}
    for ref in reference_list:
        molecule = ref['molecule']
        conjugate = ref['conjugate']
        key = (molecule, conjugate)

        if key in references:
            raise ValueError(
                f"Duplicate reference for {molecule}-{conjugate}. "
                f"Each reference set should have only one entry per molecule-conjugate pair."
            )
        references[key] = ref
    return references


def load_and_process_reference(
    directory: str,
    molecule: str,
    conjugate: str,
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
        conjugate: What the molecule is conjugated to (e.g., "EpCAM", "BSA", "IgG")
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
                'molecule': str,
                'conjugate': str
            }
    """
    # Expand tilde in directory path
    directory = os.path.expanduser(directory)

    print(f"\nProcessing: {molecule}-{conjugate}")
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

    # Add molecule and conjugate tags
    averaged['molecule'] = molecule
    averaged['conjugate'] = conjugate

    # Generate plot - flat structure with prefix
    plot_path = f"{output_dir}/reference_{molecule}_{conjugate}.png"
    print(f"Creating plot: {plot_path}")
    plot_reference(averaged, molecule=molecule, output_path=plot_path)
    print("✓ Plot saved")

    return averaged


def load_and_process_sample(
    directory: str,
    name: str,
    molecule_conjugates: list[tuple[str, str]],
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
        molecule_conjugates: List of (molecule, conjugate) tuples present in sample
                           (e.g., [("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")])
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
                'molecule_conjugates': list[tuple[str, str]],
                'name': str,
                'replicates': list[dict]  # List of individual corrected spectra
            }
    """
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
    averaged['molecule_conjugates'] = molecule_conjugates
    averaged['name'] = name

    # Add replicates to the result
    averaged['replicates'] = corrected_spectra

    # Generate plot - flat structure with prefix
    safe_name = name.replace(" ", "_").replace("/", "-")
    plot_path = f"{output_dir}/sample_{safe_name}.png"
    print(f"Creating plot: {plot_path}")

    # Extract just molecules for plotting (plotting functions expect molecule names only)
    molecules = [mol for mol, conj in molecule_conjugates]
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
                        'contributions': {(molecule, conjugate): percentage},
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
        for mol_conj in sample_data['molecule_conjugates']:
            if mol_conj not in references:
                molecule, conjugate = mol_conj
                raise ValueError(
                    f"Sample '{sample_data['name']}' requires {molecule}-{conjugate} "
                    f"but no reference was loaded for it. "
                    f"Available references: {list(references.keys())}"
                )

    # Store results for summary
    deconv_results = {}

    for sample_key, sample_data in samples.items():
        display_name = sample_data['name']
        sample_molecule_conjugates = sample_data['molecule_conjugates']
        replicates = sample_data['replicates']

        print(f"\n{display_name}:")
        print(f"  Processing {len(replicates)} replicates...")
        print(f"  Normalizing and deconvolving (range: {wavenumber_range[0]}-{wavenumber_range[1]} cm⁻¹)...")

        # Filter references to only include the ones in this sample
        sample_references = {
            mol_conj: references[mol_conj]
            for mol_conj in sample_molecule_conjugates
        }

        # Process each replicate
        replicate_results = []

        for i, replicate in enumerate(replicates):
            # Normalize this replicate
            normalized = normalize_spectra_l2(
                sample=replicate,
                references=sample_references,
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
            references=sample_references,
            wavenumber_range=wavenumber_range
        )

        # Extract just molecules for plotting (plotting functions expect molecule names only)
        molecules_only = [mol for mol, conj in sample_molecule_conjugates]

        plot_normalization(
            sample_name=f"{display_name} (averaged)",
            sample_spectrum=normalized_avg['sample'],
            reference_spectra=normalized_avg['references'],
            molecules=molecules_only,
            wavenumber_range=wavenumber_range,
            output_path=f"{output_dir}/normalization_averaged_{sample_key}.png"
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
            output_path=f"{output_dir}/deconvolution_averaged_{sample_key}.png"
        )

        # Also create original-scale deconvolution plot
        plot_deconvolution_original_scale(
            sample_name=f"{display_name} (averaged)",
            sample_spectrum=normalized_avg['sample'],
            result=deconv_avg,
            wavenumber_range=wavenumber_range,
            output_path=f"{output_dir}/deconvolution_averaged_original_scale_{sample_key}.png"
        )

        # Calculate and print summary statistics across replicates
        import numpy as np

        contributions_by_mol_conj = {mol_conj: [] for mol_conj in sample_molecule_conjugates}
        r_squared_values = []

        for result in replicate_results:
            for mol_conj in sample_molecule_conjugates:
                contributions_by_mol_conj[mol_conj].append(result['contributions'][mol_conj])
            r_squared_values.append(result['metrics']['r_squared'])

        # Print mean ± std for each molecule-conjugate
        contributions_str = ", ".join([
            f"{mol}-{conj}={np.mean(contributions_by_mol_conj[mol_conj]):.1f}±{np.std(contributions_by_mol_conj[mol_conj]):.1f}%"
            for mol_conj in sample_molecule_conjugates
            for mol, conj in [mol_conj]  # Unpack tuple for formatting
        ])
        r2_mean = np.mean(r_squared_values)
        r2_std = np.std(r_squared_values)

        print(f"  ✓ Processed {len(replicates)} replicates")
        print(f"  ✓ Contributions: {contributions_str}")
        print(f"  ✓ R²: {r2_mean:.3f}±{r2_std:.3f}")

    # Create box plots showing distribution across all samples
    print(f"\nCreating box plots...")
    boxplot_path = f"{output_dir}/deconvolution_boxplots.png"
    plot_deconvolution_boxplots(
        samples=samples,
        deconv_results=deconv_results,
        output_path=boxplot_path
    )
    print(f"✓ Box plots saved: {boxplot_path}")

    return deconv_results


def extract_peak_intensities(
    molecule: str,
    deconv_results: dict
) -> dict[str, list[float]]:
    """
    Extract peak intensities for a molecule from deconvolution results.

    For each replicate in the deconvolution results, extracts the intensity at the
    molecule's characteristic peak from the deconvolved individual contributions.
    Returns intensities grouped by conjugate type.

    Args:
        molecule: Molecule name (e.g., "MBA", "DTNB", "TFMBA")
        deconv_results: Dictionary of replicate deconvolution results (can be filtered subset)

    Returns:
        Dict mapping conjugate names to lists of peak intensities
        (e.g., {"EpCAM": [120.5, 115.3, ...], "BSA": [...]})
    """
    import numpy as np
    from .molecules import get_peak

    peak_wn = get_peak(molecule)
    conjugate_intensities = {}

    # Process all samples in the provided results
    for sample_key, replicate_results in deconv_results.items():
        # Extract peak intensity from each replicate
        for result in replicate_results:
            # Find the molecule-conjugate pair in this result's contributions
            mol_conj = None
            for key in result['contributions'].keys():
                mol, conj = key
                if mol == molecule:
                    mol_conj = key
                    conjugate = conj
                    break

            if mol_conj is None:
                continue  # This molecule not in this result

            # Initialize list for this conjugate if needed
            if conjugate not in conjugate_intensities:
                conjugate_intensities[conjugate] = []

            # Get the individual contribution (normalized scale)
            contribution_normalized = np.array(result['individual_contributions'][mol_conj])

            # Scale back to original using norm_factor
            norm_factor = result['norm_factor']
            contribution_original = contribution_normalized * norm_factor

            # Find the index closest to the peak wavenumber
            wavenumbers = np.array(result['wavenumbers'])
            peak_idx = np.argmin(np.abs(wavenumbers - peak_wn))

            # Get intensity at peak
            peak_intensity = contribution_original[peak_idx]
            conjugate_intensities[conjugate].append(peak_intensity)

    return conjugate_intensities


def calculate_histogram_scales(
    group_intensities: dict[str, dict[str, dict[str, list[float]]]],
    bin_size: int = 25
) -> tuple[float, int]:
    """
    Calculate unified x_max and y_max for histogram datasets across multiple groups.

    Computes the maximum intensity value (for x-axis) and maximum bin count (for y-axis)
    across all molecules and groups, enabling consistent scales across multiple histograms.

    Args:
        group_intensities: Nested dict structure:
            {group_prefix: {molecule: {conjugate: [intensities]}}}
        bin_size: Bin width for histogram (default 25)

    Returns:
        Tuple of (x_max, y_max) where:
            x_max: Maximum intensity value across all data
            y_max: Maximum bin count (with 10% padding) across all histograms
    """
    import numpy as np

    # Collect all intensities to find x_max
    all_intensities = []
    for group_data in group_intensities.values():
        for mol_data in group_data.values():
            for intensities in mol_data.values():
                all_intensities.extend(intensities)

    if not all_intensities:
        return (0.0, 0)

    x_max = max(all_intensities)

    # Calculate y_max by binning all datasets
    bins = np.arange(0, x_max + bin_size, bin_size)
    max_count = 0

    for group_data in group_intensities.values():
        for mol_data in group_data.values():
            for intensities in mol_data.values():
                if intensities:
                    counts, _ = np.histogram(intensities, bins=bins)
                    max_count = max(max_count, max(counts))

    y_max = int(max_count * 1.1)  # Add 10% padding

    return (x_max, y_max)


def plot_all_peak_histograms(
    deconv_results: dict,
    groups: dict[str, list[str]],
    output_dir: str,
    bin_size: int = 25
) -> None:
    """
    Plot peak intensity histograms for all molecules across multiple sample groups.

    Creates histograms showing the distribution of peak intensities for each molecule
    (MBA, DTNB, TFMBA), with unified scales across all groups for comparison.

    Args:
        deconv_results: Dictionary of deconvolution results from normalize_and_deconvolve_samples
        groups: Dictionary mapping group prefix to list of sample keys
            e.g., {"mcf7_rep1": ["MCF7_AB_1", "MCF7_BSA_1", "MCF7_IgG_1"]}
            Use empty string "" for no prefix
        output_dir: Directory to save histogram plots
        bin_size: Bin width for histograms (default 25)

    Output files:
        peak_intensity_histogram_{prefix}_{molecule}.png for each group/molecule
        (or peak_intensity_histogram_{molecule}.png if prefix is empty)
    """
    from .plotting import plot_peak_intensity_histogram
    from .molecules import get_all_molecules

    molecules = get_all_molecules()

    # 1. Extract intensities for all groups and molecules
    group_intensities = {}
    for prefix, keys in groups.items():
        selected = {k: deconv_results[k] for k in keys}
        group_intensities[prefix] = {
            mol: extract_peak_intensities(mol, selected)
            for mol in molecules
        }

    # 2. Calculate unified scales across all data
    x_max, y_max = calculate_histogram_scales(group_intensities, bin_size=bin_size)

    # 3. Plot each group's histograms with unified scales
    for prefix, mol_intensities in group_intensities.items():
        for molecule, conj_intensities in mol_intensities.items():
            # Build output filename
            if prefix:
                filename = f"peak_intensity_histogram_{prefix}_{molecule}.png"
            else:
                filename = f"peak_intensity_histogram_{molecule}.png"

            output_path = f"{output_dir}/{filename}"

            plot_peak_intensity_histogram(
                molecule=molecule,
                conjugate_intensities=conj_intensities,
                output_path=output_path,
                bin_size=bin_size,
                x_max=x_max,
                y_max=y_max
            )
            print(f"  ✓ Saved: {filename}")


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
    for (molecule, conjugate), data in sorted(references.items()):
        print(f"  {molecule}-{conjugate}: {data['count']} spectra")

    print(f"\nSamples processed:")
    for sample_key, data in sorted(samples.items(), key=lambda x: x[1]['name']):
        print(f"  {data['name']}: {data['count']} spectra")

    # Group samples by their conjugate types (e.g., all Ab samples, all BSA samples)
    import numpy as np
    from collections import defaultdict

    # Group samples by their molecule_conjugates tuple (frozen for hashing)
    groups = defaultdict(list)
    for sample_key, sample_data in samples.items():
        # Use frozenset of conjugates to group (e.g., all samples with EpCAM/HER2/TROP2)
        conjugate_set = frozenset(conj for mol, conj in sample_data['molecule_conjugates'])
        groups[conjugate_set].append(sample_key)

    # Print a table for each conjugate group
    for conjugate_set in sorted(groups.keys(), key=lambda x: sorted(x)):
        sample_keys = groups[conjugate_set]

        # Get the molecule-conjugates from first sample (all in group have same ones)
        first_sample = samples[sample_keys[0]]
        group_mol_conj = sorted(first_sample['molecule_conjugates'])

        # Determine group name from conjugates
        conjugate_list = sorted(conjugate_set)
        if len(conjugate_list) == 1:
            group_name = conjugate_list[0]
        else:
            group_name = "/".join(conjugate_list)

        print(f"\n{'='*60}")
        print(f"Deconvolution Results - {group_name} Conjugates")
        print(f"{'='*60}")

        # Build header
        header = f"{'Sample':<25}"
        for mol, conj in group_mol_conj:
            label = f"{mol}-{conj}"
            header += f" {label:>20}"
        header += f" {'R²':>14}"
        print(f"\n{header}")
        print("-" * (25 + len(group_mol_conj) * 21 + 15))

        # Print each sample in this group (sorted alphabetically by name)
        for sample_key in sorted(sample_keys, key=lambda k: samples[k]['name']):
            display_name = samples[sample_key]['name']
            replicate_results = deconv_results[sample_key]
            sample_mol_conj = samples[sample_key]['molecule_conjugates']

            # Calculate mean ± std for each molecule-conjugate
            contributions_by_mol_conj = {mol_conj: [] for mol_conj in sample_mol_conj}
            r_squared_values = []

            for result in replicate_results:
                for mol_conj in sample_mol_conj:
                    contributions_by_mol_conj[mol_conj].append(result['contributions'][mol_conj])
                r_squared_values.append(result['metrics']['r_squared'])

            # Build row
            row = f"{display_name:<25}"
            for mol_conj in group_mol_conj:
                mean = np.mean(contributions_by_mol_conj[mol_conj])
                std = np.std(contributions_by_mol_conj[mol_conj])
                row += f" {f'{mean:.1f}±{std:.1f}%':>20}"

            r2_mean = np.mean(r_squared_values)
            r2_std = np.std(r_squared_values)
            row += f" {f'{r2_mean:.3f}±{r2_std:.3f}':>14}"
            print(row)

    print(f"\nPlots saved in:")
    print(f"  {output_dir}/")
