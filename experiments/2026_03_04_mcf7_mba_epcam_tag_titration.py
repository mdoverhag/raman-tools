#!/usr/bin/env python3
"""
Experiment: MCF7 MBA EpCAM — Tag Titration (10uL, 25uL, 50uL, 65uL)
Date: 2026-03-04
Description: Singleplex MBA peak intensity analysis of fixed MCF7 cells
             with EpCAM conjugate at varying tag volumes (10, 25, 50, 65 uL).
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
    "2026-03-04 MCF7 MBA EpCAM - Tag Titration 10ul 25ul 50ul 60ul"
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
    # 10uL
    "MBA_10uL_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 10uL 1",
        name="MBA EpCAM 10uL 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10uL_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 10uL 2",
        name="MBA EpCAM 10uL 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10uL_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 10uL 3",
        name="MBA EpCAM 10uL 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 25uL
    "MBA_25uL_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 25uL 1",
        name="MBA EpCAM 25uL 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_25uL_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 25uL 2",
        name="MBA EpCAM 25uL 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_25uL_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 25uL 3",
        name="MBA EpCAM 25uL 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 50uL
    "MBA_50uL_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 50uL 1",
        name="MBA EpCAM 50uL 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_50uL_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 50uL 2",
        name="MBA EpCAM 50uL 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_50uL_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 50uL 3",
        name="MBA EpCAM 50uL 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 65uL
    "MBA_65uL_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 65uL 1",
        name="MBA EpCAM 65uL 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_65uL_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 65uL 2",
        name="MBA EpCAM 65uL 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_65uL_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7 MBA EpCAM 65uL 3",
        name="MBA EpCAM 65uL 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
}


# ============================================================
# Step 2: Load tags-only sample
# ============================================================

print("\n" + "=" * 60)
print("LOADING TAGS-ONLY SAMPLE")
print("=" * 60)

tag_samples = {
    "Tag_EpCAM": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM tags only",
        name="MBA EpCAM tags",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
}


# ============================================================
# Step 3: Peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

VOLUME_COLORS = {
    "MBA EpCAM 10uL 1": "#1f77b4", "MBA EpCAM 10uL 2": "#1f77b4", "MBA EpCAM 10uL 3": "#1f77b4",
    "MBA EpCAM 25uL 1": "#ff7f0e", "MBA EpCAM 25uL 2": "#ff7f0e", "MBA EpCAM 25uL 3": "#ff7f0e",
    "MBA EpCAM 50uL 1": "#2ca02c", "MBA EpCAM 50uL 2": "#2ca02c", "MBA EpCAM 50uL 3": "#2ca02c",
    "MBA EpCAM 65uL 1": "#d62728", "MBA EpCAM 65uL 2": "#d62728", "MBA EpCAM 65uL 3": "#d62728",
}

plot_peak_histograms_from_samples(
    mcf7_samples,
    groups={"cells_only": list(mcf7_samples.keys())},
    molecules=["MBA"],
    output_dir=output,
    title_prefix="MCF7 EpCAM Titration",
    colors=VOLUME_COLORS,
)


# ============================================================
# Summary
# ============================================================

experiment_summary(
    experiment,
    samples={**mcf7_samples, **tag_samples},
    output_dir=output,
)
