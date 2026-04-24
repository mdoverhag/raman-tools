#!/usr/bin/env python3
"""
Experiment: MCF7 and SKBR3 Specificity — MBA
Date: 2026-04-01
Description: Singleplex MBA peak intensity analysis of MCF7 and SKBR3 cells
             across four antibody conjugates (EpCAM, HER2, TROP2) with BSA
             and IgG controls. Separate histograms per cell line, unified
             scales for direct comparison.
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
    "2026-04-01-Specifity MCF7 and SKBR3"
)

# Experiment name (derived from script filename)
experiment = os.path.splitext(os.path.basename(__file__))[0]

# Create output directory (auto-versioned) in results/
output = create_output_dir(experiment, base_dir="results")
print(f"Output directory: {output}\n")


# ============================================================
# Step 1: Load samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    # MCF7
    "MCF7_EpCAM": load_and_process_sample(
        f"{DATA_DIR}/MCF7_EpCAM.txt",
        name="MCF7 EpCAM",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF7_HER2": load_and_process_sample(
        f"{DATA_DIR}/MCF7_HER2.txt",
        name="MCF7 HER2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MCF7_TROP2": load_and_process_sample(
        f"{DATA_DIR}/MCF7_TROP2.txt",
        name="MCF7 TROP2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MCF7_BSA": load_and_process_sample(
        f"{DATA_DIR}/MCF7_BSA.txt",
        name="MCF7 BSA",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MCF7_IgG": load_and_process_sample(
        f"{DATA_DIR}/MCF7_IgG.txt",
        name="MCF7 IgG",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    # SKBR3
    "SKBR3_EpCAM": load_and_process_sample(
        f"{DATA_DIR}/SKBR3_EpCAM.txt",
        name="SKBR3 EpCAM",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "SKBR3_HER2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3_HER2.txt",
        name="SKBR3 HER2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "SKBR3_TROP2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3_TROP2.txt",
        name="SKBR3 TROP2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "SKBR3_BSA": load_and_process_sample(
        f"{DATA_DIR}/SKBR3_BSA.txt",
        name="SKBR3 BSA",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "SKBR3_IgG": load_and_process_sample(
        f"{DATA_DIR}/SKBR3_IgG.txt",
        name="SKBR3 IgG",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
}


# ============================================================
# Step 2: Peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_samples(
    samples,
    groups={
        "mcf7": [
            "MCF7_EpCAM", "MCF7_HER2", "MCF7_TROP2",
            "MCF7_BSA", "MCF7_IgG",
        ],
        "skbr3": [
            "SKBR3_EpCAM", "SKBR3_HER2", "SKBR3_TROP2",
            "SKBR3_BSA", "SKBR3_IgG",
        ],
    },
    molecules=["MBA"],
    output_dir=output,
)


# ============================================================
# Summary
# ============================================================

experiment_summary(
    experiment,
    samples=samples,
    output_dir=output,
)
