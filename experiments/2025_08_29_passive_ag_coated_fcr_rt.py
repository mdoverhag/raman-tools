#!/usr/bin/env python3
"""
Experiment: MCF7 Passive Ag-coated with FcR Block (RT)
Date: 2025-08-29
Description: MCF7 passive conjugation, silver coated, with FcR Block at room temperature for 45 min.
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

# Data directories
REFERENCE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive"
)
SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-08-29 Passive Ag coated Fc receptor RT 45min"
)

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-08-29-passive-ag-coated-fcr-rt", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Step 1: Load and process references (from 2025-07-17)
# ============================================================

print("=" * 60)
print("LOADING REFERENCES")
print("=" * 60)

references = build_reference_dict(
    [
        # Antibody-conjugated references
        load_and_process_reference(
            f"{REFERENCE_DIR}/MBA EpCAM",
            molecule="MBA",
            conjugate="EpCAM",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/DTNB HER2",
            molecule="DTNB",
            conjugate="HER2",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/TFMBA TROP2",
            molecule="TFMBA",
            conjugate="TROP2",
            output_dir=output,
        ),
        # BSA control references
        load_and_process_reference(
            f"{REFERENCE_DIR}/MBA BSA",
            molecule="MBA",
            conjugate="BSA",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/DTNB BSA",
            molecule="DTNB",
            conjugate="BSA",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/TFMBA BSA",
            molecule="TFMBA",
            conjugate="BSA",
            output_dir=output,
        ),
        # IgG control references
        load_and_process_reference(
            f"{REFERENCE_DIR}/MBA IgG",
            molecule="MBA",
            conjugate="IgG",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/DTNB IgG",
            molecule="DTNB",
            conjugate="IgG",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/TFMBA IgG",
            molecule="TFMBA",
            conjugate="IgG",
            output_dir=output,
        ),
    ]
)


# ============================================================
# Step 2: Load and process multiplex samples
# ============================================================

print("\n" + "=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    "Multi_Ab_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi Ab 1",
        name="Multi Ab 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Multi_Ab_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi Ab 2",
        name="Multi Ab 2",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Multi_Ab_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi Ab 3",
        name="Multi Ab 3",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Multi_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi BSA 1",
        name="Multi BSA 1",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Multi_BSA_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi BSA 2",
        name="Multi BSA 2",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Multi_BSA_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi BSA 3",
        name="Multi BSA 3",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Multi_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi IgG 1",
        name="Multi IgG 1",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "Multi_IgG_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi IgG 2",
        name="Multi IgG 2",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "Multi_IgG_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Multi IgG 3",
        name="Multi IgG 3",
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
    groups={
        "rep1": ["Multi_Ab_1", "Multi_BSA_1", "Multi_IgG_1"],
        "rep2": ["Multi_Ab_2", "Multi_BSA_2", "Multi_IgG_2"],
        "rep3": ["Multi_Ab_3", "Multi_BSA_3", "Multi_IgG_3"],
    },
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
