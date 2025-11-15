#!/usr/bin/env python3
"""
Experiment: SKBR3 Passive Analysis
Date: 2025-07-17
Description: Process reference spectra for DTNB/HER2, MBA/EpCAM, and TFMBA/TROP2
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples
)

# Data directory
DATA_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive")

# Wavenumber range for normalization
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-07-17-skbr3-passive", base_dir="results")
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

print("="*60)
print("LOADING SAMPLES")
print("="*60)

samples = {
    "Multiplex_Ab": load_and_process_sample(
        f"{DATA_DIR}/Multiplex Ab",
        name="Multiplex Ab",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
}

print(f"✓ Processed {len(samples)} sample\n")


# ============================================================
# Step 3 & 4: Normalize and deconvolute each sample
# ============================================================

deconv_results = normalize_and_deconvolve_samples(
    samples=samples,
    references=references,
    wavenumber_range=WAVENUMBER_RANGE,
    output_dir=output,
    molecules=["MBA", "DTNB", "TFMBA"]
)


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
for sample_key, data in samples.items():
    print(f"  {data['name']}: {data['count']} spectra")

print(f"\nDeconvolution Results:")
print(f"\n{'Sample':<25} {'MBA':>8} {'DTNB':>8} {'TFMBA':>8} {'R²':>8}")
print("-" * 60)
for sample_key, result in deconv_results.items():
    display_name = samples[sample_key]['name']
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
