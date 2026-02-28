#!/usr/bin/env python3
"""
Experiment: MCF7 MBA Fixed Cells — EpCAM
Date: 2026-02-19
Description: Singleplex MBA peak intensity analysis of fixed MCF7 cells
             with EpCAM antibody conjugate, BSA and IgG controls.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
    plot_peak_histograms_from_samples,
    experiment_summary,
)

# Data directory
DATA_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/"
    "2026-02-19 MCF7 MBA fixed cells"
)

# Experiment name (derived from script filename)
experiment = os.path.splitext(os.path.basename(__file__))[0]

# Create output directory (auto-versioned) in results/
output = create_output_dir(experiment, base_dir="results")
print(f"Output directory: {output}\n")


# ============================================================
# Step 1: Load MCF7 cell samples
# ============================================================

print("=" * 60)
print("LOADING MCF7 SAMPLES")
print("=" * 60)

mcf7_samples = {
    # EpCAM
    "MBA_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/EpCAM_1.txt",
        name="MBA EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/EpCAM_2.txt",
        name="MBA EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/EpCAM_3.txt",
        name="MBA EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # BSA controls
    "MBA_BSA_1": load_and_process_sample(
        f"{DATA_DIR}/BSA_1.txt",
        name="MBA BSA 1",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_BSA_2": load_and_process_sample(
        f"{DATA_DIR}/BSA_2.txt",
        name="MBA BSA 2",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_BSA_3": load_and_process_sample(
        f"{DATA_DIR}/BSA_3.txt",
        name="MBA BSA 3",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    # IgG controls
    "MBA_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/IgG_1.txt",
        name="MBA IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MBA_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/IgG_2.txt",
        name="MBA IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MBA_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/IgG_3.txt",
        name="MBA IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
}


# ============================================================
# Step 2: Load tags-only samples
# ============================================================

print("\n" + "=" * 60)
print("LOADING TAGS-ONLY SAMPLES")
print("=" * 60)

tag_samples = {
    "Tag_EpCAM": load_and_process_sample(
        f"{DATA_DIR}/EpCAM_tags only.txt",
        name="MBA EpCAM tags",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Tag_BSA": load_and_process_sample(
        f"{DATA_DIR}/BSA_tags only.txt",
        name="MBA BSA tags",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Tag_IgG": load_and_process_sample(
        f"{DATA_DIR}/IgG_tags only.txt",
        name="MBA IgG tags",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
}


# ============================================================
# Step 3: MCF7 peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("MCF7 PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_samples(
    mcf7_samples,
    groups={
        "epcam": [
            "MBA_EpCAM_1", "MBA_EpCAM_2", "MBA_EpCAM_3",
            "MBA_IgG_1", "MBA_IgG_2", "MBA_IgG_3",
            "MBA_BSA_1", "MBA_BSA_2", "MBA_BSA_3",
        ],
    },
    molecules=["MBA"],
    output_dir=output,
)


# ============================================================
# Step 4: Tags-only peak intensity histogram
# ============================================================

print("\n" + "=" * 60)
print("TAGS-ONLY PEAK INTENSITY HISTOGRAM")
print("=" * 60)

plot_peak_histograms_from_samples(
    tag_samples,
    groups={
        "tags": ["Tag_EpCAM", "Tag_BSA", "Tag_IgG"],
    },
    molecules=["MBA"],
    output_dir=output,
    bin_size=250,
)


# ============================================================
# Summary
# ============================================================

experiment_summary(
    experiment,
    samples={**mcf7_samples, **tag_samples},
    output_dir=output,
)
