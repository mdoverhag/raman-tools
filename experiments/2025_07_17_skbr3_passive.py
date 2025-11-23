#!/usr/bin/env python3
"""
Experiment: SKBR3 Passive Analysis
Date: 2025-07-17
Description: Process reference spectra for DTNB/HER2, MBA/EpCAM, and TFMBA/TROP2
"""

import os
from raman_lib import (
    create_output_dir,
    build_reference_dict,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples,
    extract_peak_intensities,
    plot_peak_intensity_histogram,
    print_experiment_summary,
)

# Data directory
DATA_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive"
)

# Wavenumber range for normalization
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-07-17-skbr3-passive", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Step 1: Load and process references
# ============================================================

print("=" * 60)
print("LOADING REFERENCES")
print("=" * 60)

references = build_reference_dict([
    # Antibody-conjugated references
    load_and_process_reference(
        f"{DATA_DIR}/MBA EpCAM", molecule="MBA", conjugate="EpCAM", output_dir=output
    ),
    load_and_process_reference(
        f"{DATA_DIR}/DTNB HER2", molecule="DTNB", conjugate="HER2", output_dir=output
    ),
    load_and_process_reference(
        f"{DATA_DIR}/TFMBA TROP2", molecule="TFMBA", conjugate="TROP2", output_dir=output
    ),
    # BSA control references
    load_and_process_reference(
        f"{DATA_DIR}/MBA BSA", molecule="MBA", conjugate="BSA", output_dir=output
    ),
    load_and_process_reference(
        f"{DATA_DIR}/DTNB BSA", molecule="DTNB", conjugate="BSA", output_dir=output
    ),
    load_and_process_reference(
        f"{DATA_DIR}/TFMBA BSA", molecule="TFMBA", conjugate="BSA", output_dir=output
    ),
    # IgG control references
    load_and_process_reference(
        f"{DATA_DIR}/MBA IgG", molecule="MBA", conjugate="IgG", output_dir=output
    ),
    load_and_process_reference(
        f"{DATA_DIR}/DTNB IgG", molecule="DTNB", conjugate="IgG", output_dir=output
    ),
    load_and_process_reference(
        f"{DATA_DIR}/TFMBA IgG", molecule="TFMBA", conjugate="IgG", output_dir=output
    ),
])


# ============================================================
# Step 2: Load and process samples
# ============================================================

print("\n" + "=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    "Multiplex_Ab": load_and_process_sample(
        f"{DATA_DIR}/Multiplex Ab",
        name="Multiplex Ab",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Multiplex_BSA": load_and_process_sample(
        f"{DATA_DIR}/Multiplex BSA",
        name="Multiplex BSA",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Multiplex_IgG": load_and_process_sample(
        f"{DATA_DIR}/Multiplex IgG",
        name="Multiplex IgG",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
}


# ============================================================
# Step 3: Normalize and deconvolve all replicates
# ============================================================

print("\n" + "=" * 60)
print("NORMALIZATION & DECONVOLUTION")
print("=" * 60)

deconv_results = normalize_and_deconvolve_samples(
    samples=samples,
    references=references,
    wavenumber_range=WAVENUMBER_RANGE,
    output_dir=output,
)


# ============================================================
# Peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

# Filter to just the samples we want to plot
selected_results = {
    k: deconv_results[k]
    for k in ["Multiplex_Ab", "Multiplex_BSA", "Multiplex_IgG"]
}

# Extract intensities for all three molecules
mba_intensities = extract_peak_intensities(
    molecule="MBA",
    deconv_results=selected_results
)
dtnb_intensities = extract_peak_intensities(
    molecule="DTNB",
    deconv_results=selected_results
)
tfmba_intensities = extract_peak_intensities(
    molecule="TFMBA",
    deconv_results=selected_results
)

# Calculate global max intensity (x-axis)
all_intensities = []
for dataset in [mba_intensities, dtnb_intensities, tfmba_intensities]:
    for intensities in dataset.values():
        all_intensities.extend(intensities)
x_max = max(all_intensities)

# Calculate global max count (y-axis) by binning all datasets
import numpy as np
bin_size = 25
bins = np.arange(0, x_max + bin_size, bin_size)
max_count = 0

for dataset in [mba_intensities, dtnb_intensities, tfmba_intensities]:
    for intensities in dataset.values():
        counts, _ = np.histogram(intensities, bins=bins)
        max_count = max(max_count, max(counts))

y_max = int(max_count * 1.1)  # Add 10% padding

# Plot all three with same scales
plot_peak_intensity_histogram(
    molecule="MBA",
    conjugate_intensities=mba_intensities,
    output_path=f"{output}/peak_intensity_histogram_MBA.png",
    x_max=x_max,
    y_max=y_max
)
print("  ✓ Saved: peak_intensity_histogram_MBA.png")

plot_peak_intensity_histogram(
    molecule="DTNB",
    conjugate_intensities=dtnb_intensities,
    output_path=f"{output}/peak_intensity_histogram_DTNB.png",
    x_max=x_max,
    y_max=y_max
)
print("  ✓ Saved: peak_intensity_histogram_DTNB.png")

plot_peak_intensity_histogram(
    molecule="TFMBA",
    conjugate_intensities=tfmba_intensities,
    output_path=f"{output}/peak_intensity_histogram_TFMBA.png",
    x_max=x_max,
    y_max=y_max
)
print("  ✓ Saved: peak_intensity_histogram_TFMBA.png")


# ============================================================
# Summary
# ============================================================

print("\n" + "=" * 60)
print("EXPERIMENT COMPLETE")
print("=" * 60)

print_experiment_summary(
    output_dir=output,
    references=references,
    samples=samples,
    deconv_results=deconv_results,
)
