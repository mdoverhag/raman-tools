#!/usr/bin/env python3
"""
Experiment: Linker Bare vs Ag-coated AuNP EpCAM
Date: 2025-08-01
Description: Linker Antibody conjugation Bare vs Ag-coated AuNP with MCF7 cells.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
    plot_peak_histograms_from_samples,
)

# Data directories
AG_SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-08-01 Linker Ag coated EpCAM only"
)
BARE_SAMPLE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-08-01 Linker Bare AuNP EpCAM only"
)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-08-01-linker-bare-vs-ag-coated", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Load and process samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    # Ag-coated samples
    "Ag_EpCAM_1": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA EpCAM 1",
        name="Ag-coated EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Ag_EpCAM_2": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA EpCAM 2",
        name="Ag-coated EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Ag_EpCAM_3": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA EpCAM 3",
        name="Ag-coated EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Ag_BSA_1": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA BSA 1",
        name="Ag-coated BSA 1",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Ag_BSA_2": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA BSA 2",
        name="Ag-coated BSA 2",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Ag_BSA_3": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA BSA 3",
        name="Ag-coated BSA 3",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Ag_IgG_1": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA IgG 1",
        name="Ag-coated IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "Ag_IgG_2": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA IgG 2",
        name="Ag-coated IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "Ag_IgG_3": load_and_process_sample(
        f"{AG_SAMPLE_DIR}/MBA IgG 3",
        name="Ag-coated IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    # Bare samples
    "Bare_EpCAM_1": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA EpCAM 1",
        name="Bare EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Bare_EpCAM_2": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA EpCAM 2",
        name="Bare EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Bare_EpCAM_3": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA EpCAM 3",
        name="Bare EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Bare_BSA_1": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA BSA 1",
        name="Bare BSA 1",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Bare_BSA_2": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA BSA 2",
        name="Bare BSA 2",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Bare_BSA_3": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA BSA 3",
        name="Bare BSA 3",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Bare_IgG_1": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA IgG 1",
        name="Bare IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "Bare_IgG_2": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA IgG 2",
        name="Bare IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "Bare_IgG_3": load_and_process_sample(
        f"{BARE_SAMPLE_DIR}/Bare MBA IgG 3",
        name="Bare IgG 3",
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
        "ag_rep1": ["Ag_EpCAM_1", "Ag_BSA_1", "Ag_IgG_1"],
        "ag_rep2": ["Ag_EpCAM_2", "Ag_BSA_2", "Ag_IgG_2"],
        "ag_rep3": ["Ag_EpCAM_3", "Ag_BSA_3", "Ag_IgG_3"],
        "bare_rep1": ["Bare_EpCAM_1", "Bare_BSA_1", "Bare_IgG_1"],
        "bare_rep2": ["Bare_EpCAM_2", "Bare_BSA_2", "Bare_IgG_2"],
        "bare_rep3": ["Bare_EpCAM_3", "Bare_BSA_3", "Bare_IgG_3"],
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
