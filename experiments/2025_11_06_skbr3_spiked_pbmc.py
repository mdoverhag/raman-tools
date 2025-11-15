#!/usr/bin/env python3
"""
Experiment: SKBR3 Spiked PBMCs Multiplex Analysis
Date: 2025-11-06
Description: Process three multiplex Ab samples using references from 2025-07-17
"""

import os
from raman_lib import (
    normalize_spectra_l2,
    plot_normalization,
    deconvolve_nnls,
    plot_deconvolution,
    create_output_dir,
    ensure_output_subdir,
    load_and_process_reference,
    load_and_process_sample
)

# Data directories
REFERENCE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive")
SAMPLE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-11-06 SKBR3 spiked PBMCs multiplex")

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-11-06-skbr3-spiked-pbmc", base_dir="results")
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
# Step 2: Load and process multiplex samples (from 2025-11-06)
# ============================================================

print("="*60)
print("LOADING SAMPLES")
print("="*60)

samples = {
    "Multiplex_Ab_1": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 1",
        name="Multiplex Ab 1",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "Multiplex_Ab_2": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 2",
        name="Multiplex Ab 2",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "Multiplex_Ab_3": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 3",
        name="Multiplex Ab 3",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
}

print(f"✓ Processed {len(samples)} samples\n")


# ============================================================
# Step 3 & 4: Normalize and deconvolute each sample
# ============================================================

# Create subdirectories for normalization and deconvolution plots
norm_dir = ensure_output_subdir(output, "normalization")
deconv_dir = ensure_output_subdir(output, "deconvolution")

for sample_key, sample_data in samples.items():
    sample_name = sample_key.replace("_", " ")

    # Step 3: Normalization
    print("="*60)
    print(f"Normalizing: {sample_name}")
    print("="*60)

    print(f"Normalization range: {WAVENUMBER_RANGE[0]}-{WAVENUMBER_RANGE[1]} cm⁻¹")
    print("Applying L2 normalization...")

    normalized = normalize_spectra_l2(
        sample=sample_data,
        references=references,
        wavenumber_range=WAVENUMBER_RANGE
    )

    print("✓ Normalization complete")

    print(f"Creating plot: {norm_dir}/{sample_key}.png")
    plot_normalization(
        sample_name=sample_name,
        sample_spectrum=normalized['sample'],
        reference_spectra=normalized['references'],
        molecules=["MBA", "DTNB", "TFMBA"],
        wavenumber_range=WAVENUMBER_RANGE,
        output_path=f"{norm_dir}/{sample_key}.png"
    )
    print("✓ Plot saved\n")

    # Step 4: Deconvolution
    print("="*60)
    print(f"Deconvoluting: {sample_name}")
    print("="*60)

    print(f"Analysis range: {WAVENUMBER_RANGE[0]}-{WAVENUMBER_RANGE[1]} cm⁻¹")
    print("Performing NNLS deconvolution...")

    deconv_result = deconvolve_nnls(
        sample_spectrum=normalized['sample'],
        reference_spectra=normalized['references'],
        wavenumber_range=WAVENUMBER_RANGE
    )

    print("✓ Deconvolution complete")

    print("\nContributions:")
    for molecule, percentage in sorted(deconv_result['contributions'].items()):
        print(f"  {molecule}: {percentage:.1f}%")

    print(f"\nMetrics:")
    print(f"  RMSE: {deconv_result['metrics']['rmse']:.3f}")
    print(f"  R²: {deconv_result['metrics']['r_squared']:.3f}")

    print(f"\nCreating plot: {deconv_dir}/{sample_key}.png")
    plot_deconvolution(
        sample_name=sample_name,
        sample_spectrum=normalized['sample'],
        result=deconv_result,
        wavenumber_range=WAVENUMBER_RANGE,
        output_path=f"{deconv_dir}/{sample_key}.png"
    )
    print("✓ Plot saved\n")


# ============================================================
# Summary
# ============================================================

print("="*60)
print("EXPERIMENT COMPLETE")
print("="*60)
print(f"Output directory: {output}")
print(f"\nProcessed references:")
for molecule, data in references.items():
    print(f"  {molecule}: {data['count']} spectra averaged")
print(f"\nProcessed samples:")
for sample_name, data in samples.items():
    print(f"  {sample_name}: {data['count']} spectra averaged")
print(f"\nPlots saved:")
print(f"  References: {output}/references/")
print(f"  Samples: {output}/samples/")
print(f"  Normalization: {output}/normalization/")
print(f"  Deconvolution: {output}/deconvolution/")
print("="*60)
