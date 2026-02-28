#!/usr/bin/env python3
"""
Experiment: MCF7 Double-coated MBA EpCAM vs BSA
Date: 2025-12-05
Description: Single-plex MBA-EpCAM analysis of MCF7 cells using double-coated AuNPs.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
    plot_peak_histograms_from_samples,
    experiment_summary,
)

# Data directory
SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-12-05 AuNP Ag MCF7 double coated MBA EpCAM vs BSA"
)

# Experiment name (derived from script filename)
experiment = os.path.splitext(os.path.basename(__file__))[0]

# Create output directory (auto-versioned) in results/
output = create_output_dir(experiment, base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Load and process samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    "MCF7_EpCAM_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag double coated MBA EpCAM",
        name="MCF7 EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF7_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag double coated MBA BSA",
        name="MCF7 BSA 1",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
}

# ============================================================
# Peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_samples(
    samples,
    groups={
        "": ["MCF7_EpCAM_1", "MCF7_BSA_1"],
    },
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
