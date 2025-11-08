#!/usr/bin/env python3
"""
Experiment: SKBR3 Spiked PBMCs Multiplex Analysis
Date: 2025-11-06
Description: Process three multiplex Ab samples using references from 2025-07-17
"""

import os
from raman_lib import (
    load_spectra,
    apply_baseline_correction,
    calculate_average,
    plot_reference,
    plot_sample,
    normalize_spectra_l2,
    plot_normalization,
    deconvolve_nnls,
    plot_deconvolution,
    create_output_dir,
    ensure_output_subdir
)

# Data directories
REFERENCE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive")
SAMPLE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-11-06 SKBR3 spiked PBMCs multiplex")

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("skbr3-spiked-pbmc-2025-11-06", base_dir="results")
print(f"Output directory: {output}\n")

# Create subdirectory for reference plots
refs_dir = ensure_output_subdir(output, "references")


# ============================================================
# Step 1: Load and process references (from 2025-07-17)
# ============================================================

references = {}

# Reference 1: DTNB HER2
print("="*60)
print("Processing: DTNB HER2 (Reference)")
print("="*60)

dtnb_dir = f"{REFERENCE_DIR}/DTNB HER2"
print(f"Loading from: {dtnb_dir}")

dtnb_spectra = load_spectra(dtnb_dir)
print(f"✓ Loaded {len(dtnb_spectra)} spectra")

print("Applying baseline correction...")
dtnb_corrected = [apply_baseline_correction(s) for s in dtnb_spectra]
print("✓ Baseline correction applied")

print("Calculating average...")
dtnb_averaged = calculate_average(dtnb_corrected)
print(f"✓ Average calculated from {dtnb_averaged['count']} spectra")

print(f"Creating plot: {refs_dir}/DTNB.png")
plot_reference(dtnb_averaged, molecule="DTNB", output_path=f"{refs_dir}/DTNB.png")
print("✓ Plot saved\n")

references["DTNB"] = dtnb_averaged


# Reference 2: MBA EpCAM
print("="*60)
print("Processing: MBA EpCAM (Reference)")
print("="*60)

mba_dir = f"{REFERENCE_DIR}/MBA EpCAM"
print(f"Loading from: {mba_dir}")

mba_spectra = load_spectra(mba_dir)
print(f"✓ Loaded {len(mba_spectra)} spectra")

print("Applying baseline correction...")
mba_corrected = [apply_baseline_correction(s) for s in mba_spectra]
print("✓ Baseline correction applied")

print("Calculating average...")
mba_averaged = calculate_average(mba_corrected)
print(f"✓ Average calculated from {mba_averaged['count']} spectra")

print(f"Creating plot: {refs_dir}/MBA.png")
plot_reference(mba_averaged, molecule="MBA", output_path=f"{refs_dir}/MBA.png")
print("✓ Plot saved\n")

references["MBA"] = mba_averaged


# Reference 3: TFMBA TROP2
print("="*60)
print("Processing: TFMBA TROP2 (Reference)")
print("="*60)

tfmba_dir = f"{REFERENCE_DIR}/TFMBA TROP2"
print(f"Loading from: {tfmba_dir}")

tfmba_spectra = load_spectra(tfmba_dir)
print(f"✓ Loaded {len(tfmba_spectra)} spectra")

print("Applying baseline correction...")
tfmba_corrected = [apply_baseline_correction(s) for s in tfmba_spectra]
print("✓ Baseline correction applied")

print("Calculating average...")
tfmba_averaged = calculate_average(tfmba_corrected)
print(f"✓ Average calculated from {tfmba_averaged['count']} spectra")

print(f"Creating plot: {refs_dir}/TFMBA.png")
plot_reference(tfmba_averaged, molecule="TFMBA", output_path=f"{refs_dir}/TFMBA.png")
print("✓ Plot saved\n")

references["TFMBA"] = tfmba_averaged


# ============================================================
# Step 2: Load and process multiplex samples (from 2025-11-06)
# ============================================================

samples = {}

# Create subdirectory for sample plots
samples_dir = ensure_output_subdir(output, "samples")

# Sample 1: Multiplex Ab 1
print("="*60)
print("Processing: Multiplex Ab 1")
print("="*60)

multiplex_ab1_dir = f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 1"
print(f"Loading from: {multiplex_ab1_dir}")

multiplex_ab1_spectra = load_spectra(multiplex_ab1_dir)
print(f"✓ Loaded {len(multiplex_ab1_spectra)} spectra")

print("Applying baseline correction...")
multiplex_ab1_corrected = [apply_baseline_correction(s) for s in multiplex_ab1_spectra]
print("✓ Baseline correction applied")

print("Calculating average...")
multiplex_ab1_averaged = calculate_average(multiplex_ab1_corrected)
print(f"✓ Average calculated from {multiplex_ab1_averaged['count']} spectra")

print(f"Creating plot: {samples_dir}/Multiplex_Ab_1.png")
plot_sample(
    multiplex_ab1_averaged,
    sample_name="Multiplex Ab 1",
    molecules=["MBA", "DTNB", "TFMBA"],
    output_path=f"{samples_dir}/Multiplex_Ab_1.png"
)
print("✓ Plot saved\n")

samples["Multiplex_Ab_1"] = multiplex_ab1_averaged


# Sample 2: Multiplex Ab 2
print("="*60)
print("Processing: Multiplex Ab 2")
print("="*60)

multiplex_ab2_dir = f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 2"
print(f"Loading from: {multiplex_ab2_dir}")

multiplex_ab2_spectra = load_spectra(multiplex_ab2_dir)
print(f"✓ Loaded {len(multiplex_ab2_spectra)} spectra")

print("Applying baseline correction...")
multiplex_ab2_corrected = [apply_baseline_correction(s) for s in multiplex_ab2_spectra]
print("✓ Baseline correction applied")

print("Calculating average...")
multiplex_ab2_averaged = calculate_average(multiplex_ab2_corrected)
print(f"✓ Average calculated from {multiplex_ab2_averaged['count']} spectra")

print(f"Creating plot: {samples_dir}/Multiplex_Ab_2.png")
plot_sample(
    multiplex_ab2_averaged,
    sample_name="Multiplex Ab 2",
    molecules=["MBA", "DTNB", "TFMBA"],
    output_path=f"{samples_dir}/Multiplex_Ab_2.png"
)
print("✓ Plot saved\n")

samples["Multiplex_Ab_2"] = multiplex_ab2_averaged


# Sample 3: Multiplex Ab 3
print("="*60)
print("Processing: Multiplex Ab 3")
print("="*60)

multiplex_ab3_dir = f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 3"
print(f"Loading from: {multiplex_ab3_dir}")

multiplex_ab3_spectra = load_spectra(multiplex_ab3_dir)
print(f"✓ Loaded {len(multiplex_ab3_spectra)} spectra")

print("Applying baseline correction...")
multiplex_ab3_corrected = [apply_baseline_correction(s) for s in multiplex_ab3_spectra]
print("✓ Baseline correction applied")

print("Calculating average...")
multiplex_ab3_averaged = calculate_average(multiplex_ab3_corrected)
print(f"✓ Average calculated from {multiplex_ab3_averaged['count']} spectra")

print(f"Creating plot: {samples_dir}/Multiplex_Ab_3.png")
plot_sample(
    multiplex_ab3_averaged,
    sample_name="Multiplex Ab 3",
    molecules=["MBA", "DTNB", "TFMBA"],
    output_path=f"{samples_dir}/Multiplex_Ab_3.png"
)
print("✓ Plot saved\n")

samples["Multiplex_Ab_3"] = multiplex_ab3_averaged


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
print(f"  References: {refs_dir}/")
print(f"  Samples: {samples_dir}/")
print(f"  Normalization: {norm_dir}/")
print(f"  Deconvolution: {deconv_dir}/")
print("="*60)
