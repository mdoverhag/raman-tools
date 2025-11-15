#!/usr/bin/env python3
"""
Experiment: MCF7 vs SKBR3 Spiked PBMCs Multiplex Analysis
Date: 2025-11-04
Description: Process all multiplex Ab samples (3x MCF7, 3x SKBR3) using references from 2025-07-17
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples
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

print("\n" + "="*60)
print("LOADING SAMPLES")
print("="*60)

samples = {
    "MCF7_Ab_1": load_and_process_sample(
        f"{SAMPLE_DIR}/MCF7 spiked PBMC_10E5_Multiplex Ab 1",
        name="MCF7 Multiplex Ab 1",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "MCF7_Ab_2": load_and_process_sample(
        f"{SAMPLE_DIR}/MCF7 spiked PBMC_10E5_Multiplex Ab 2",
        name="MCF7 Multiplex Ab 2",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "MCF7_Ab_3": load_and_process_sample(
        f"{SAMPLE_DIR}/MCF7 spiked PBMC_10E5_Multiplex Ab 3",
        name="MCF7 Multiplex Ab 3",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "SKBR3_Ab_1": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 1",
        name="SKBR3 Multiplex Ab 1",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "SKBR3_Ab_2": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 2",
        name="SKBR3 Multiplex Ab 2",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "SKBR3_Ab_3": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 3",
        name="SKBR3 Multiplex Ab 3",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
}

print(f"✓ Processed {len(samples)} samples\n")


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
