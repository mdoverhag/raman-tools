#!/usr/bin/env python3
"""
Experiment: MCF7 MBA EpCAM — LOD Cell Titrations (10, 100, 1000, 10000, 100000 cells)
Date: 2026-03-10
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
    "2026-03-10 MCF7 - MBA EpCAM - LOD"
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
    # 10 cells
    "MBA_10cells_1": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 10 cells 1",
        name="MBA EpCAM 10 cells 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10cells_2": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 10 cells 2",
        name="MBA EpCAM 10 cells 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10cells_3": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 10 cells 3",
        name="MBA EpCAM 10 cells 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 100 cells
    "MBA_100cells_1": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 100 cells 1",
        name="MBA EpCAM 100 cells 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_100cells_2": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 100 cells 2",
        name="MBA EpCAM 100 cells 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_100cells_3": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 100 cells 3",
        name="MBA EpCAM 100 cells 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 1000 cells
    "MBA_1000cells_1": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 1000 cells 1",
        name="MBA EpCAM 1000 cells 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_1000cells_2": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 1000 cells 2",
        name="MBA EpCAM 1000 cells 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_1000cells_3": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 1000 cells 3",
        name="MBA EpCAM 1000 cells 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 10000 cells
    "MBA_10000cells_1": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 10000 cells 1",
        name="MBA EpCAM 10000 cells 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10000cells_2": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 10000 cells 2",
        name="MBA EpCAM 10000 cells 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10000cells_3": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 10000 cells 3",
        name="MBA EpCAM 10000 cells 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # 100000 cells
    "MBA_100000cells_1": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 100 000 cells 1",
        name="MBA EpCAM 100000 cells 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_100000cells_2": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 100 000 cells 2",
        name="MBA EpCAM 100000 cells 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_100000cells_3": load_and_process_sample(
        f"{DATA_DIR}/LOD MCF7 - MBA EpCAM - 100 000 cells 3",
        name="MBA EpCAM 100000 cells 3",
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
    "MBA EpCAM 10 cells 1": "#1f77b4", "MBA EpCAM 10 cells 2": "#1f77b4", "MBA EpCAM 10 cells 3": "#1f77b4",
    "MBA EpCAM 100 cells 1": "#ff7f0e", "MBA EpCAM 100 cells 2": "#ff7f0e", "MBA EpCAM 100 cells 3": "#ff7f0e",
    "MBA EpCAM 1000 cells 1": "#2ca02c", "MBA EpCAM 1000 cells 2": "#2ca02c", "MBA EpCAM 1000 cells 3": "#2ca02c",
    "MBA EpCAM 10000 cells 1": "#d62728", "MBA EpCAM 10000 cells 2": "#d62728", "MBA EpCAM 10000 cells 3": "#d62728",
    "MBA EpCAM 100000 cells 1": "#9467bd", "MBA EpCAM 100000 cells 2": "#9467bd", "MBA EpCAM 100000 cells 3": "#9467bd",
}

plot_peak_histograms_from_samples(
    samples,
    groups={"all": list(samples.keys())},
    molecules=["MBA"],
    output_dir=output,
    title_prefix="MCF7 EpCAM LOD",
    colors=CONC_COLORS,
)


# ============================================================
# Summary
# ============================================================

experiment_summary(
    experiment,
    samples=samples,
    output_dir=output,
)
