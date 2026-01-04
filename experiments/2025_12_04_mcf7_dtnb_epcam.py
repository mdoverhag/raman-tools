#!/usr/bin/env python3
"""
Experiment: MCF7 DTNB-EpCAM Single-plex
Date: 2025-12-04
Description: Single-plex DTNB-EpCAM analysis of MCF7 cells.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
    plot_peak_histograms_from_samples,
)

# Data directory
SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-12-04 AuNP Ag MCF7 DTNB-EpCAM only"
)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-12-04-mcf7-dtnb-epcam", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Load and process samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    "MCF7_EpCAM_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MCF7 DTNB-EpCAM",
        name="MCF7 EpCAM 1",
        molecule_conjugates=[("DTNB", "EpCAM")],
        output_dir=output,
    ),
    "MCF7_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/DTNB BSA1",
        name="MCF7 BSA 1",
        molecule_conjugates=[("DTNB", "BSA")],
        output_dir=output,
    ),
    "MCF7_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/DTNB IgG1",
        name="MCF7 IgG 1",
        molecule_conjugates=[("DTNB", "IgG")],
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
        "": ["MCF7_EpCAM_1", "MCF7_BSA_1", "MCF7_IgG_1"],
    },
    output_dir=output,
    x_max=500,
)


# ============================================================
# Summary
# ============================================================

print("\n" + "=" * 60)
print("EXPERIMENT COMPLETE")
print("=" * 60)

print(f"\nOutput directory: {output}")
print(f"\nSamples processed:")
for sample_key, data in sorted(samples.items(), key=lambda x: x[1]["name"]):
    print(f"  {data['name']}: {data['count']} spectra")
