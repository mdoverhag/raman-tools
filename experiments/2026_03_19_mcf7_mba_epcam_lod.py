#!/usr/bin/env python3
"""
Experiment: MCF7 MBA EpCAM — Limit of Detection
Date: 2026-03-19
Description: Singleplex MBA peak intensity analysis of MCF7 cells with EpCAM
             conjugate at varying cell concentrations (10 to 100000 cells) to
             determine limit of detection (LOD).
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
    "2026-03-19 MCF7 - EpCAM - Limit of Detection"
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
    "MBA_10cells": load_and_process_sample(
        f"{DATA_DIR}/LOD_1.txt",
        name="MBA EpCAM 10 cells",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_100cells": load_and_process_sample(
        f"{DATA_DIR}/LOD_2.txt",
        name="MBA EpCAM 100 cells",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_1000cells": load_and_process_sample(
        f"{DATA_DIR}/LOD_3.txt",
        name="MBA EpCAM 1000 cells",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_10000cells": load_and_process_sample(
        f"{DATA_DIR}/LOD_4.txt",
        name="MBA EpCAM 10000 cells",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_100000cells": load_and_process_sample(
        f"{DATA_DIR}/LOD_5.txt",
        name="MBA EpCAM 100000 cells",
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
    "MBA EpCAM 10 cells": "#1f77b4",
    "MBA EpCAM 100 cells": "#ff7f0e",
    "MBA EpCAM 1000 cells": "#2ca02c",
    "MBA EpCAM 10000 cells": "#d62728",
    "MBA EpCAM 100000 cells": "#9467bd",
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
