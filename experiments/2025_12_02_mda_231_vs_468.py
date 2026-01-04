#!/usr/bin/env python3
"""
Experiment: MDA-MB-231 vs MDA-MB-468 PD-L1 Analysis
Date: 2025-12-02
Description: Single-plex TFMBA-PD-L1 analysis comparing MDA-MB-231 and MDA-MB-468 cell lines.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
    plot_peak_histograms_from_samples,
)

# Data directory
SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-12-02 AuNP Ag MDA MB 231 vs 468 TFMBA PD-L1"
)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-12-02-mda-231-vs-468", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Load and process samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    # MDA-MB-231 samples
    "MDA231_PDL1_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA MB 231 PD-L1 1",
        name="MDA-231 PD-L1 1",
        molecule_conjugates=[("TFMBA", "PD-L1")],
        output_dir=output,
    ),
    "MDA231_PDL1_2": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA MB 231 PD-L1 2",
        name="MDA-231 PD-L1 2",
        molecule_conjugates=[("TFMBA", "PD-L1")],
        output_dir=output,
    ),
    "MDA231_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 231 BSA1",
        name="MDA-231 BSA 1",
        molecule_conjugates=[("TFMBA", "BSA")],
        output_dir=output,
    ),
    "MDA231_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 231 IgG1",
        name="MDA-231 IgG 1",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
    ),
    # MDA-MB-468 samples
    "MDA468_PDL1_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 468 PD-L1 1",
        name="MDA-468 PD-L1 1",
        molecule_conjugates=[("TFMBA", "PD-L1")],
        output_dir=output,
    ),
    "MDA468_BSA_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 468 BSA1",
        name="MDA-468 BSA 1",
        molecule_conjugates=[("TFMBA", "BSA")],
        output_dir=output,
    ),
    "MDA468_IgG_1": load_and_process_sample(
        f"{SAMPLE_DIR}/AuNP Ag MDA 468 IgG1",
        name="MDA-468 IgG 1",
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
        "mda231": ["MDA231_PDL1_1", "MDA231_PDL1_2", "MDA231_BSA_1", "MDA231_IgG_1"],
        "mda468": ["MDA468_PDL1_1", "MDA468_BSA_1", "MDA468_IgG_1"],
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
for sample_key, data in sorted(samples.items(), key=lambda x: x[1]['name']):
    print(f"  {data['name']}: {data['count']} spectra")
