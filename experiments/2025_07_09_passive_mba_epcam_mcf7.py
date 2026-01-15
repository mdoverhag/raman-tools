#!/usr/bin/env python3
"""
Experiment: MCF7 Passive MBA-EpCAM Single-plex
Date: 2025-07-09
Description: Single-plex MBA-EpCAM analysis of MCF7 cells with passive antibody conjugation.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
    plot_peak_histograms_from_samples,
)

# Data directory
SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-07-09 passive Ab conjugation, MBA EpCAM, MCF7 with Controls"
)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-07-09-passive-mba-epcam-mcf7", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Load and process samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    "MBA_EpCAM_1": load_and_process_sample(
        f"{SAMPLE_DIR}/MBA EpCAM 1",
        name="MBA EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_EpCAM_2": load_and_process_sample(
        f"{SAMPLE_DIR}/MBA EpCAM 2",
        name="MBA EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/MBA BSA 1",
        name="MBA BSA 1",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_BSA_2": load_and_process_sample(
        f"{SAMPLE_DIR}/MBA BSA 2",
        name="MBA BSA 2",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/MBA IgG 1",
        name="MBA IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MBA_IgG_2": load_and_process_sample(
        f"{SAMPLE_DIR}/MBA IgG 2",
        name="MBA IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
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
        "rep1": ["MBA_EpCAM_1", "MBA_BSA_1", "MBA_IgG_1"],
        "rep2": ["MBA_EpCAM_2", "MBA_BSA_2", "MBA_IgG_2"],
    },
    output_dir=output,
    x_max=2000,
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
