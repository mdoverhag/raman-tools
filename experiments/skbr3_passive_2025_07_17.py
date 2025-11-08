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
    plot_reference,
    plot_sample,
    create_output_dir,
    ensure_output_subdir
)

# Data directory
DATA_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive")

# Create output directory (auto-versioned) in results/
output = create_output_dir("skbr3-passive-2025-07-17", base_dir="results")
print(f"Output directory: {output}\n")

# Create subdirectory for reference plots
refs_dir = ensure_output_subdir(output, "references")

# ============================================================
# Step 1: Load and process references
# ============================================================

references = {}

# Reference 1: DTNB HER2
print("="*60)
print("Processing: DTNB HER2")
print("="*60)

dtnb_dir = f"{DATA_DIR}/DTNB HER2"
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
print("Processing: MBA EpCAM")
print("="*60)

mba_dir = f"{DATA_DIR}/MBA EpCAM"
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
print("Processing: TFMBA TROP2")
print("="*60)

tfmba_dir = f"{DATA_DIR}/TFMBA TROP2"
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
print("="*60)
