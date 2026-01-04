#!/usr/bin/env python3
"""
Experiment: MCF7 vs SKBR3 Multiplex Analysis
Date: 2025-11-25
Description: MCF7 and SKBR3 spiked in healthy donor PBMCs. 
"""

import os
from raman_lib import (
    create_output_dir,
    build_reference_dict,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples,
    print_experiment_summary,
)

# Data directories
REFERENCE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive"
)
SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-11-25 AuNP Ag MCF7 vs SKBR3 repeat Ep H2 T2"
)

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-11-25-skbr3-vs-mcf7", base_dir="results")
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
# Step 2: Load and process multiplex samples (from 2025-11-25)
# ============================================================

print("\n" + "=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    "MCF7_AB_1": load_and_process_sample(
        f"{SAMPLE_DIR}/MCF7 Ab1",
        name="MCF7 AB 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "MCF7_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/MCF7 BSA1",
        name="MCF7 BSA 1",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "MCF7_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/MCF7 IgG1",
        name="MCF7 IgG 1",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "SKBR3_AB_2": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 Ab 2",
        name="SKBR3 AB 2",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "SKBR3_BSA_2": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 BSA 2",
        name="SKBR3 BSA 2",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "SKBR3_IgG_2": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 IgG2",
        name="SKBR3 IgG 2",
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
