#!/usr/bin/env python3
"""
Experiment: MDA-MB-231 and MDA-MB-468 Specificity - MBA
Date: 2026-06-04
Description: Singleplex MBA peak intensity analysis of MDA-MB-231 and
             MDA-MB-468 cells across EpCAM, HER2, PD-L1, and TROP2 antibody
             conjugates with IgG controls. Histograms compare each target's
             three replicates against matched IgG replicates for each cell line.
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
    "2026-06-04 MDA-MB-231 MDA-MB-468 - Ep H2 PD T2 IgG"
)

# Experiment name (derived from script filename)
experiment = os.path.splitext(os.path.basename(__file__))[0]

# Create output directory (auto-versioned) in results/
output = create_output_dir(experiment, base_dir="results")
print(f"Output directory: {output}\n")

RAW_OUTPUT_DIR = os.path.join(output, "loaded_samples")


# ============================================================
# Step 1: Load samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    # MDA-MB-231
    "MDA231_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-EpCAM 1",
        name="MDA-MB-231 EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-EpCAM 2",
        name="MDA-MB-231 EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-EpCAM 3",
        name="MDA-MB-231 EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_HER2_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-HER2 1",
        name="MDA-MB-231 HER2 1",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_HER2_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-HER2 2",
        name="MDA-MB-231 HER2 2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_HER2_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-HER2 3",
        name="MDA-MB-231 HER2 3",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-PDL1 1",
        name="MDA-MB-231 PD-L1 1",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-PDL1 2",
        name="MDA-MB-231 PD-L1 2",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_PDL1_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-PDL1 3",
        name="MDA-MB-231 PD-L1 3",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_TROP2_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-TROP2 1",
        name="MDA-MB-231 TROP2 1",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_TROP2_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-TROP2 2",
        name="MDA-MB-231 TROP2 2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_TROP2_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-TROP2 3",
        name="MDA-MB-231 TROP2 3",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-IgG 1",
        name="MDA-MB-231 IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-IgG 2",
        name="MDA-MB-231 IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231/MDA-MB-231 - MBA-IgG 3",
        name="MDA-MB-231 IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    # MDA-MB-468
    "MDA468_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-EpCAM 1",
        name="MDA-MB-468 EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-EpCAM 2",
        name="MDA-MB-468 EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-EpCAM 3",
        name="MDA-MB-468 EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_HER2_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-HER2 1",
        name="MDA-MB-468 HER2 1",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_HER2_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-HER2 2",
        name="MDA-MB-468 HER2 2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_HER2_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-HER2 3",
        name="MDA-MB-468 HER2 3",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-PDL1 1",
        name="MDA-MB-468 PD-L1 1",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-PDL1 2",
        name="MDA-MB-468 PD-L1 2",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_PDL1_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-PDL1 3",
        name="MDA-MB-468 PD-L1 3",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_TROP2_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-TROP2 1",
        name="MDA-MB-468 TROP2 1",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_TROP2_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-TROP2 2",
        name="MDA-MB-468 TROP2 2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_TROP2_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-TROP2 3",
        name="MDA-MB-468 TROP2 3",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-IgG 1",
        name="MDA-MB-468 IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-IgG 2",
        name="MDA-MB-468 IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468/MDA-MB-468 - MBA-IgG 3",
        name="MDA-MB-468 IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
}


# ============================================================
# Step 2: Peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_samples(
    samples,
    groups={
        "mda231_epcam": [
            "MDA231_EpCAM_1", "MDA231_EpCAM_2", "MDA231_EpCAM_3",
            "MDA231_IgG_1", "MDA231_IgG_2", "MDA231_IgG_3",
        ],
        "mda231_her2": [
            "MDA231_HER2_1", "MDA231_HER2_2", "MDA231_HER2_3",
            "MDA231_IgG_1", "MDA231_IgG_2", "MDA231_IgG_3",
        ],
        "mda231_pdl1": [
            "MDA231_PDL1_1", "MDA231_PDL1_2", "MDA231_PDL1_3",
            "MDA231_IgG_1", "MDA231_IgG_2", "MDA231_IgG_3",
        ],
        "mda231_trop2": [
            "MDA231_TROP2_1", "MDA231_TROP2_2", "MDA231_TROP2_3",
            "MDA231_IgG_1", "MDA231_IgG_2", "MDA231_IgG_3",
        ],
        "mda468_epcam": [
            "MDA468_EpCAM_1", "MDA468_EpCAM_2", "MDA468_EpCAM_3",
            "MDA468_IgG_1", "MDA468_IgG_2", "MDA468_IgG_3",
        ],
        "mda468_her2": [
            "MDA468_HER2_1", "MDA468_HER2_2", "MDA468_HER2_3",
            "MDA468_IgG_1", "MDA468_IgG_2", "MDA468_IgG_3",
        ],
        "mda468_pdl1": [
            "MDA468_PDL1_1", "MDA468_PDL1_2", "MDA468_PDL1_3",
            "MDA468_IgG_1", "MDA468_IgG_2", "MDA468_IgG_3",
        ],
        "mda468_trop2": [
            "MDA468_TROP2_1", "MDA468_TROP2_2", "MDA468_TROP2_3",
            "MDA468_IgG_1", "MDA468_IgG_2", "MDA468_IgG_3",
        ],
    },
    molecules=["MBA"],
    x_max=325,
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
