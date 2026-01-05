#!/usr/bin/env python3
"""
Experiment: MCF7 Passive Analysis
Date: 2025-07-21
Description: MCF7 Passive Analysis
"""

import os
from raman_lib import (
    create_output_dir,
    build_reference_dict,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples,
    plot_peak_histograms_from_deconv,
    print_experiment_summary,
)

# Data directory
DATA_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-07-21 MCF7 passive"
)

# Wavenumber range for normalization
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-07-21-mcf7-passive", base_dir="results")
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
    "Multi_Ab": load_and_process_sample(
        f"{DATA_DIR}/Multi Ab",
        name="Multi Ab",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Multi_BSA": load_and_process_sample(
        f"{DATA_DIR}/Multi BSA",
        name="Multi BSA",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Multi_IgG": load_and_process_sample(
        f"{DATA_DIR}/Multi IgG",
        name="Multi IgG",
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

plot_peak_histograms_from_deconv(
    deconv_results,
    groups={"": ["Multi_Ab", "Multi_BSA", "Multi_IgG"]},
    output_dir=output,
)


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
