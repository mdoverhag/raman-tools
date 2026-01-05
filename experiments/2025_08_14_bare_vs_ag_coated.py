#!/usr/bin/env python3
"""
Experiment: MCF7 Bare vs Ag-coated AuNP Multiplex Analysis
Date: 2025-08-14
Description: Comparing bare vs silver-coated gold nanoparticles on MCF7 cells with multiplex antibodies.
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
    "~/Documents/Spectroscopy Results/2025-08-14 Passive Bare vs Ag coated AuNP"
)

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-08-14-bare-vs-ag-coated", base_dir="results")
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
    # Ag-coated samples
    "Ag_AB_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi Ab 1",
        name="Ag-coated AB 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Ag_AB_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi Ab 2",
        name="Ag-coated AB 2",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Ag_AB_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi Ab 3",
        name="Ag-coated AB 3",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Ag_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi BSA 1",
        name="Ag-coated BSA 1",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Ag_BSA_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi BSA 2",
        name="Ag-coated BSA 2",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Ag_BSA_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi BSA 3",
        name="Ag-coated BSA 3",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Ag_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi IgG 1",
        name="Ag-coated IgG 1",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "Ag_IgG_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi IgG 2",
        name="Ag-coated IgG 2",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "Ag_IgG_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Ag coated Multi IgG 3",
        name="Ag-coated IgG 3",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    # Bare samples
    "Bare_AB_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi Ab 1",
        name="Bare AB 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Bare_AB_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi Ab 2",
        name="Bare AB 2",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Bare_AB_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi Ab 3",
        name="Bare AB 3",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "Bare_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi BSA 1",
        name="Bare BSA 1",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Bare_BSA_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi BSA 2",
        name="Bare BSA 2",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Bare_BSA_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi BSA 3",
        name="Bare BSA 3",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "Bare_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi IgG 1",
        name="Bare IgG 1",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "Bare_IgG_2": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi IgG 2",
        name="Bare IgG 2",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "Bare_IgG_3": load_and_process_sample(
        f"{SAMPLE_DIR}/Bare Multi IgG 3",
        name="Bare IgG 3",
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
        "ag_rep1": ["Ag_AB_1", "Ag_BSA_1", "Ag_IgG_1"],
        "ag_rep2": ["Ag_AB_2", "Ag_BSA_2", "Ag_IgG_2"],
        "ag_rep3": ["Ag_AB_3", "Ag_BSA_3", "Ag_IgG_3"],
        "bare_rep1": ["Bare_AB_1", "Bare_BSA_1", "Bare_IgG_1"],
        "bare_rep2": ["Bare_AB_2", "Bare_BSA_2", "Bare_IgG_2"],
        "bare_rep3": ["Bare_AB_3", "Bare_BSA_3", "Bare_IgG_3"],
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
