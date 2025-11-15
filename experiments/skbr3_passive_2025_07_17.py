#!/usr/bin/env python3
"""
Experiment: SKBR3 Passive Analysis
Date: 2025-07-17
Description: Process reference spectra for DTNB/HER2, MBA/EpCAM, and TFMBA/TROP2
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

# Data directory
DATA_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive")

# Wavenumber range for normalization
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("skbr3-passive-2025-07-17", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Step 1: Load and process references
# ============================================================

print("="*60)
print("LOADING REFERENCES")
print("="*60)

references = {
    "MBA": load_and_process_reference(
        f"{DATA_DIR}/MBA EpCAM",
        molecule="MBA",
        output_dir=output
    ),
    "DTNB": load_and_process_reference(
        f"{DATA_DIR}/DTNB HER2",
        molecule="DTNB",
        output_dir=output
    ),
    "TFMBA": load_and_process_reference(
        f"{DATA_DIR}/TFMBA TROP2",
        molecule="TFMBA",
        output_dir=output
    ),
}

print(f"✓ Processed {len(references)} references\n")


# ============================================================
# Step 2: Load and process samples
# ============================================================

samples = {}

# Create subdirectory for sample plots
samples_dir = ensure_output_subdir(output, "samples")

# Sample 1: Multiplex Ab
print("="*60)
print("Processing: Multiplex Ab")
print("="*60)

multiplex_ab_dir = f"{DATA_DIR}/Multiplex Ab"
print(f"Loading from: {multiplex_ab_dir}")

multiplex_ab_spectra = load_spectra(multiplex_ab_dir)
print(f"✓ Loaded {len(multiplex_ab_spectra)} spectra")

print("Applying baseline correction...")
multiplex_ab_corrected = [apply_baseline_correction(s) for s in multiplex_ab_spectra]
print("✓ Baseline correction applied")

print("Calculating average...")
multiplex_ab_averaged = calculate_average(multiplex_ab_corrected)
print(f"✓ Average calculated from {multiplex_ab_averaged['count']} spectra")

print(f"Creating plot: {samples_dir}/Multiplex_Ab.png")
plot_sample(
    multiplex_ab_averaged,
    sample_name="Multiplex Ab",
    molecules=["MBA", "DTNB", "TFMBA"],
    output_path=f"{samples_dir}/Multiplex_Ab.png"
)
print("✓ Plot saved\n")

samples["Multiplex_Ab"] = multiplex_ab_averaged


# ============================================================
# Step 3: Normalize sample against references
# ============================================================

# Create subdirectory for normalization plots
norm_dir = ensure_output_subdir(output, "normalization")

print("="*60)
print("Normalizing: Multiplex Ab")
print("="*60)

print(f"Normalization range: {WAVENUMBER_RANGE[0]}-{WAVENUMBER_RANGE[1]} cm⁻¹")
print("Applying L2 normalization...")

# Normalize sample and references
normalized = normalize_spectra_l2(
    sample=multiplex_ab_averaged,
    references=references,
    wavenumber_range=WAVENUMBER_RANGE
)

print("✓ Normalization complete")

# Create normalization plot
print(f"Creating plot: {norm_dir}/Multiplex_Ab.png")
plot_normalization(
    sample_name="Multiplex Ab",
    sample_spectrum=normalized['sample'],
    reference_spectra=normalized['references'],
    molecules=["MBA", "DTNB", "TFMBA"],
    wavenumber_range=WAVENUMBER_RANGE,
    output_path=f"{norm_dir}/Multiplex_Ab.png"
)
print("✓ Plot saved\n")


# ============================================================
# Step 4: Deconvolute multiplex spectrum
# ============================================================

# Create subdirectory for deconvolution plots
deconv_dir = ensure_output_subdir(output, "deconvolution")

print("="*60)
print("Deconvoluting: Multiplex Ab")
print("="*60)

print(f"Analysis range: {WAVENUMBER_RANGE[0]}-{WAVENUMBER_RANGE[1]} cm⁻¹")
print("Performing NNLS deconvolution...")

# Deconvolute the normalized sample against normalized references
deconv_result = deconvolve_nnls(
    sample_spectrum=normalized['sample'],
    reference_spectra=normalized['references'],
    wavenumber_range=WAVENUMBER_RANGE
)

print("✓ Deconvolution complete")

# Display results
print("\nContributions:")
for molecule, percentage in sorted(deconv_result['contributions'].items()):
    print(f"  {molecule}: {percentage:.1f}%")

print(f"\nMetrics:")
print(f"  RMSE: {deconv_result['metrics']['rmse']:.3f}")
print(f"  R²: {deconv_result['metrics']['r_squared']:.3f}")

# Create deconvolution plot
print(f"\nCreating plot: {deconv_dir}/Multiplex_Ab.png")
plot_deconvolution(
    sample_name="Multiplex Ab",
    sample_spectrum=normalized['sample'],
    result=deconv_result,
    wavenumber_range=WAVENUMBER_RANGE,
    output_path=f"{deconv_dir}/Multiplex_Ab.png"
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
