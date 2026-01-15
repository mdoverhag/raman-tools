#!/usr/bin/env python3
"""
Experiment: RT vs 37°C Incubation Temperature Comparison
Date: 2025-09-23
Description: Comparing multiplex antibody binding at room temperature (21°C) vs 37°C incubation.
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
    "~/Documents/Spectroscopy Results/2025-09-23 RT vs 37degrees"
)

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-09-23-rt-vs-37degrees", base_dir="results")
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
    # 37°C samples
    "37deg_Ab_1": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi Ab 1",
        name="37°C Ab 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "37deg_Ab_2": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi Ab 2",
        name="37°C Ab 2",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "37deg_Ab_3": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi Ab 3",
        name="37°C Ab 3",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "37deg_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi BSA 1",
        name="37°C BSA 1",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "37deg_BSA_2": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi BSA 2",
        name="37°C BSA 2",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "37deg_BSA_3": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi BSA 3",
        name="37°C BSA 3",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "37deg_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi IgG 1",
        name="37°C IgG 1",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "37deg_IgG_2": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi IgG 2",
        name="37°C IgG 2",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "37deg_IgG_3": load_and_process_sample(
        f"{SAMPLE_DIR}/37deg Multi IgG 3",
        name="37°C IgG 3",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    # RT (21°C) samples
    "RT_Ab_1": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi Ab 1",
        name="RT Ab 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "RT_Ab_2": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi Ab 2",
        name="RT Ab 2",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "RT_Ab_3": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi Ab 3",
        name="RT Ab 3",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output,
    ),
    "RT_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi BSA 1",
        name="RT BSA 1",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "RT_BSA_2": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi BSA 2",
        name="RT BSA 2",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "RT_BSA_3": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi BSA 3",
        name="RT BSA 3",
        molecule_conjugates=[("MBA", "BSA"), ("DTNB", "BSA"), ("TFMBA", "BSA")],
        output_dir=output,
    ),
    "RT_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi IgG 1",
        name="RT IgG 1",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "RT_IgG_2": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi IgG 2",
        name="RT IgG 2",
        molecule_conjugates=[("MBA", "IgG"), ("DTNB", "IgG"), ("TFMBA", "IgG")],
        output_dir=output,
    ),
    "RT_IgG_3": load_and_process_sample(
        f"{SAMPLE_DIR}/RT21Deg Multi IgG 3",
        name="RT IgG 3",
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
        "37deg_rep1": ["37deg_Ab_1", "37deg_BSA_1", "37deg_IgG_1"],
        "37deg_rep2": ["37deg_Ab_2", "37deg_BSA_2", "37deg_IgG_2"],
        "37deg_rep3": ["37deg_Ab_3", "37deg_BSA_3", "37deg_IgG_3"],
        "rt_rep1": ["RT_Ab_1", "RT_BSA_1", "RT_IgG_1"],
        "rt_rep2": ["RT_Ab_2", "RT_BSA_2", "RT_IgG_2"],
        "rt_rep3": ["RT_Ab_3", "RT_BSA_3", "RT_IgG_3"],
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
