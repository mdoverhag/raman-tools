#!/usr/bin/env python3
"""
Experiment: MDA-MB-231 MBA — PD-L1 and EpCAM
Date: 2026-03-19
Description: Singleplex MBA peak intensity analysis of MDA-MB-231 cells
             with PD-L1 and EpCAM antibody conjugates, with IgG controls.
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
    "2026-03-19 MDA MB 231 - PD-L1 and EpCAM"
)

# Experiment name (derived from script filename)
experiment = os.path.splitext(os.path.basename(__file__))[0]

# Create output directory (auto-versioned) in results/
output = create_output_dir(experiment, base_dir="results")
print(f"Output directory: {output}\n")


# ============================================================
# Step 1: Load MDA-MB-231 cell samples
# ============================================================

print("=" * 60)
print("LOADING MDA-MB-231 SAMPLES")
print("=" * 60)

samples = {
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
    # PD-L1
    "MBA_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/PDL1_1.txt",
        name="MBA PD-L1 1",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MBA_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/PDL1_2.txt",
        name="MBA PD-L1 2",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MBA_PDL1_3": load_and_process_sample(
        f"{DATA_DIR}/PDL1_3.txt",
        name="MBA PD-L1 3",
        molecule_conjugates=[("MBA", "PD-L1")],
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
# Step 2: Peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_samples(
    samples,
    groups={
        "all": [
            "MBA_EpCAM_1", "MBA_EpCAM_2", "MBA_EpCAM_3",
            "MBA_PDL1_1", "MBA_PDL1_2", "MBA_PDL1_3",
            "MBA_IgG_1", "MBA_IgG_2", "MBA_IgG_3",
        ],
    },
    molecules=["MBA"],
    output_dir=output,
    title_prefix="MDA-MB-231",
    x_max=500,
)


# ============================================================
# Summary
# ============================================================

experiment_summary(
    experiment,
    samples=samples,
    output_dir=output,
)
