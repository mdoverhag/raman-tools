#!/usr/bin/env python3
"""
Experiment: MCF7 vs SKBR3 Spiked PBMCs Multiplex Analysis
Date: 2025-11-04
Description: Process all multiplex Ab samples (3x MCF7, 3x SKBR3) using references from 2025-07-17
"""

import os
from raman_lib import (
    load_spectra,
    apply_baseline_correction,
    calculate_average,
    plot_sample,
    normalize_spectra_l2,
    plot_normalization,
    deconvolve_nnls,
    plot_deconvolution,
    create_output_dir,
    ensure_output_subdir,
    load_and_process_reference
)

# Data directories
REFERENCE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive")
SAMPLE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-11-04 MCF7 vs SKBR3 spiked PBMCs Multiplex")

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-11-04-mcf7-vs-skbr3", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Step 1: Load and process references (from 2025-07-17)
# ============================================================

print("="*60)
print("LOADING REFERENCES")
print("="*60)

references = {
    "MBA": load_and_process_reference(
        f"{REFERENCE_DIR}/MBA EpCAM",
        molecule="MBA",
        output_dir=output
    ),
    "DTNB": load_and_process_reference(
        f"{REFERENCE_DIR}/DTNB HER2",
        molecule="DTNB",
        output_dir=output
    ),
    "TFMBA": load_and_process_reference(
        f"{REFERENCE_DIR}/TFMBA TROP2",
        molecule="TFMBA",
        output_dir=output
    ),
}

print(f"✓ Processed {len(references)} references\n")


# ============================================================
# Step 2: Load and process multiplex samples
# ============================================================

# Define all multiplex Ab samples to process
sample_definitions = [
    ("MCF7_Ab_1", "MCF7 spiked PBMC_10E5_Multiplex Ab 1", "MCF7 Multiplex Ab 1"),
    ("MCF7_Ab_2", "MCF7 spiked PBMC_10E5_Multiplex Ab 2", "MCF7 Multiplex Ab 2"),
    ("MCF7_Ab_3", "MCF7 spiked PBMC_10E5_Multiplex Ab 3", "MCF7 Multiplex Ab 3"),
    ("SKBR3_Ab_1", "SKBR3 spiked PBMC_10E5_Multiplex Ab 1", "SKBR3 Multiplex Ab 1"),
    ("SKBR3_Ab_2", "SKBR3 spiked PBMC_10E5_Multiplex Ab 2", "SKBR3 Multiplex Ab 2"),
    ("SKBR3_Ab_3", "SKBR3 spiked PBMC_10E5_Multiplex Ab 3", "SKBR3 Multiplex Ab 3"),
]

samples = {}
samples_dir = ensure_output_subdir(output, "samples")

print("\n" + "="*60)
print("LOADING SAMPLES")
print("="*60)

for sample_key, dir_name, display_name in sample_definitions:
    print(f"\nProcessing: {display_name}")
    sample_path = f"{SAMPLE_DIR}/{dir_name}"

    spectra = load_spectra(sample_path)
    print(f"✓ Loaded {len(spectra)} spectra")

    corrected = [apply_baseline_correction(s) for s in spectra]
    averaged = calculate_average(corrected)

    plot_sample(
        averaged,
        sample_name=display_name,
        molecules=["MBA", "DTNB", "TFMBA"],
        output_path=f"{samples_dir}/{sample_key}.png"
    )
    print(f"✓ Processed and saved")

    samples[sample_key] = averaged


# ============================================================
# Step 3 & 4: Normalize and deconvolute each sample
# ============================================================

norm_dir = ensure_output_subdir(output, "normalization")
deconv_dir = ensure_output_subdir(output, "deconvolution")

print("\n" + "="*60)
print("NORMALIZATION & DECONVOLUTION")
print("="*60)

# Store results for summary
deconv_results = {}

for sample_key, _, display_name in sample_definitions:
    sample_data = samples[sample_key]

    print(f"\n{display_name}:")
    print(f"  Normalizing (range: {WAVENUMBER_RANGE[0]}-{WAVENUMBER_RANGE[1]} cm⁻¹)...")

    # Normalize
    normalized = normalize_spectra_l2(
        sample=sample_data,
        references=references,
        wavenumber_range=WAVENUMBER_RANGE
    )

    plot_normalization(
        sample_name=display_name,
        sample_spectrum=normalized['sample'],
        reference_spectra=normalized['references'],
        molecules=["MBA", "DTNB", "TFMBA"],
        wavenumber_range=WAVENUMBER_RANGE,
        output_path=f"{norm_dir}/{sample_key}.png"
    )
    print(f"  ✓ Normalized and saved")

    # Deconvolute
    print(f"  Deconvoluting...")
    deconv_result = deconvolve_nnls(
        sample_spectrum=normalized['sample'],
        reference_spectra=normalized['references'],
        wavenumber_range=WAVENUMBER_RANGE
    )

    plot_deconvolution(
        sample_name=display_name,
        sample_spectrum=normalized['sample'],
        result=deconv_result,
        wavenumber_range=WAVENUMBER_RANGE,
        output_path=f"{deconv_dir}/{sample_key}.png"
    )

    # Store results
    deconv_results[sample_key] = deconv_result

    # Print contributions
    print(f"  ✓ Deconvoluted: MBA={deconv_result['contributions']['MBA']:.1f}%, " +
          f"DTNB={deconv_result['contributions']['DTNB']:.1f}%, " +
          f"TFMBA={deconv_result['contributions']['TFMBA']:.1f}% " +
          f"(R²={deconv_result['metrics']['r_squared']:.3f})")


# ============================================================
# Summary
# ============================================================

print("\n" + "="*60)
print("EXPERIMENT COMPLETE")
print("="*60)
print(f"\nOutput directory: {output}")

print(f"\nReferences processed:")
for molecule, data in references.items():
    print(f"  {molecule}: {data['count']} spectra")

print(f"\nSamples processed:")
for sample_key, _, display_name in sample_definitions:
    data = samples[sample_key]
    print(f"  {display_name}: {data['count']} spectra")

print(f"\nDeconvolution Results:")
print(f"\n{'Sample':<25} {'MBA':>8} {'DTNB':>8} {'TFMBA':>8} {'R²':>8}")
print("-" * 60)
for sample_key, _, display_name in sample_definitions:
    result = deconv_results[sample_key]
    print(f"{display_name:<25} {result['contributions']['MBA']:>7.1f}% " +
          f"{result['contributions']['DTNB']:>7.1f}% " +
          f"{result['contributions']['TFMBA']:>7.1f}% " +
          f"{result['metrics']['r_squared']:>8.3f}")

print(f"\nPlots saved in:")
print(f"  {output}/references/")
print(f"  {output}/samples/")
print(f"  {output}/normalization/")
print(f"  {output}/deconvolution/")
print("="*60)
