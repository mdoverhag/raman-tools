#!/usr/bin/env python3
"""
Experiment: MCF7 MBA EpCAM — LOD Cell Titrations (10e3, 10e4, 10e5, 10e6)
Date: 2026-03-05
Description: Singleplex MBA peak intensity analysis of MCF7 cells with EpCAM
             conjugate at varying cell concentrations to determine limit of
             detection (LOD).
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
    "2026-03-05 MCF7 MBA EpCAM - LOD - Cell Titrations"
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

samples = {
    # 10e3
    "MBA_10e3_1": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e3 1",
        name="MBA EpCAM 10e3 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e3_2": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e3 2",
        name="MBA EpCAM 10e3 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e3_3": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e3 3",
        name="MBA EpCAM 10e3 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 10e4
    "MBA_10e4_1": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e4 1",
        name="MBA EpCAM 10e4 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e4_2": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e4 2",
        name="MBA EpCAM 10e4 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e4_3": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e4 3",
        name="MBA EpCAM 10e4 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 10e5
    "MBA_10e5_1": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e5 1",
        name="MBA EpCAM 10e5 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e5_2": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e5 2",
        name="MBA EpCAM 10e5 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e5_3": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e5 3",
        name="MBA EpCAM 10e5 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 10e6
    "MBA_10e6_1": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e6 1",
        name="MBA EpCAM 10e6 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e6_2": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e6 2",
        name="MBA EpCAM 10e6 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10e6_3": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM 10e6 3",
        name="MBA EpCAM 10e6 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
}


# ============================================================
# Step 2: Peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

CONC_COLORS = {
    "MBA EpCAM 10e3 1": "#1f77b4", "MBA EpCAM 10e3 2": "#1f77b4", "MBA EpCAM 10e3 3": "#1f77b4",
    "MBA EpCAM 10e4 1": "#ff7f0e", "MBA EpCAM 10e4 2": "#ff7f0e", "MBA EpCAM 10e4 3": "#ff7f0e",
    "MBA EpCAM 10e5 1": "#2ca02c", "MBA EpCAM 10e5 2": "#2ca02c", "MBA EpCAM 10e5 3": "#2ca02c",
    "MBA EpCAM 10e6 1": "#d62728", "MBA EpCAM 10e6 2": "#d62728", "MBA EpCAM 10e6 3": "#d62728",
}

plot_peak_histograms_from_samples(
    samples,
    groups={"all": list(samples.keys())},
    molecules=["MBA"],
    output_dir=output,
    title_prefix="MCF7 EpCAM LOD",
    colors=CONC_COLORS,
    bin_size=250,
)


# ============================================================
# Summary
# ============================================================

experiment_summary(
    experiment,
    samples=samples,
    output_dir=output,
)
