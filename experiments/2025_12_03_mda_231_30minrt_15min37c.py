#!/usr/bin/env python3
"""
Experiment: MDA-MB-231 PD-L1 30min RT + 15min 37°C
Date: 2025-12-03
Description: Single-plex TFMBA-PD-L1 analysis of MDA-MB-231 cells with 30min RT + 15min 37°C
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
    "~/Documents/Spectroscopy Results/2025-12-03 MDA 231 30minRT 15min37degrees"
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
    "MDA231_PDL1_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 231 PD-L1 1 30minRT 15min37",
        name="MDA-231 PD-L1 1",
        molecule_conjugates=[("TFMBA", "PD-L1")],
        output_dir=output,
    ),
    "MDA231_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 231 BSA 1 30minRT 15min37",
        name="MDA-231 BSA 1",
        molecule_conjugates=[("TFMBA", "BSA")],
        output_dir=output,
    ),
    "MDA231_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 231 IgG 1 30minRT 15min37",
        name="MDA-231 IgG 1",
        molecule_conjugates=[("TFMBA", "IgG")],
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
        "": ["MDA231_PDL1_1", "MDA231_BSA_1", "MDA231_IgG_1"],
    },
    output_dir=output,
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
