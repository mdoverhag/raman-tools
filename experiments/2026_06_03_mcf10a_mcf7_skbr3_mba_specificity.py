#!/usr/bin/env python3
"""
Experiment: MCF10A, MCF7, and SKBR3 Specificity - MBA
Date: 2026-06-03
Description: Singleplex MBA peak intensity analysis of MCF10A, MCF7, and
             SKBR3 cells across EpCAM, HER2, PD-L1, and TROP2 antibody
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
    "2026-06-03 MCF10A MCF7 SKBR3 - Ep H2 PD T2 IgG"
)

# Experiment name (derived from script filename)
experiment = os.path.splitext(os.path.basename(__file__))[0]

# Create output directory (auto-versioned) in results/
output = create_output_dir(experiment, base_dir="results")
print(f"Output directory: {output}\n")


# ============================================================
# Step 1: Load samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    # MCF10A
    "MCF10A_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-EpCAM 1",
        name="MCF10A EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF10A_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-EpCAM 2",
        name="MCF10A EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF10A_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-EpCAM 3",
        name="MCF10A EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF10A_HER2_1": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-HER2 1",
        name="MCF10A HER2 1",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MCF10A_HER2_2": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-HER2 2",
        name="MCF10A HER2 2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MCF10A_HER2_3": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-HER2 3",
        name="MCF10A HER2 3",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MCF10A_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-PDL1 1",
        name="MCF10A PD-L1 1",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MCF10A_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-PDL1 2",
        name="MCF10A PD-L1 2",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MCF10A_PDL1_3": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-PDL1 3",
        name="MCF10A PD-L1 3",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MCF10A_TROP2_1": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-TROP2 1",
        name="MCF10A TROP2 1",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MCF10A_TROP2_2": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-TROP2 2",
        name="MCF10A TROP2 2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MCF10A_TROP2_3": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-TROP2 3",
        name="MCF10A TROP2 3",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MCF10A_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-IgG 1",
        name="MCF10A IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MCF10A_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-IgG 2",
        name="MCF10A IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MCF10A_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/MCF10A/MCF10A - MBA-IgG 3",
        name="MCF10A IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    # MCF7
    "MCF7_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-EpCAM 1",
        name="MCF7 EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF7_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-EpCAM 2",
        name="MCF7 EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF7_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-EpCAM 3",
        name="MCF7 EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MCF7_HER2_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-HER2 1",
        name="MCF7 HER2 1",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MCF7_HER2_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-HER2 2",
        name="MCF7 HER2 2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MCF7_HER2_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-HER2 3",
        name="MCF7 HER2 3",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MCF7_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-PDL1 1",
        name="MCF7 PD-L1 1",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MCF7_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-PDL1 2",
        name="MCF7 PD-L1 2",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MCF7_PDL1_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-PDL1 3",
        name="MCF7 PD-L1 3",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "MCF7_TROP2_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-TROP2 1",
        name="MCF7 TROP2 1",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MCF7_TROP2_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-TROP2 2",
        name="MCF7 TROP2 2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MCF7_TROP2_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-TROP2 3",
        name="MCF7 TROP2 3",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MCF7_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-IgG 1",
        name="MCF7 IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MCF7_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-IgG 2",
        name="MCF7 IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MCF7_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/MCF7/MCF7 - MBA-IgG 3",
        name="MCF7 IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    # SKBR3
    "SKBR3_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-EpCAM 1",
        name="SKBR3 EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "SKBR3_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-EpCAM 2",
        name="SKBR3 EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "SKBR3_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-EpCAM 3",
        name="SKBR3 EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "SKBR3_HER2_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-HER2 1",
        name="SKBR3 HER2 1",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "SKBR3_HER2_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-HER2 2",
        name="SKBR3 HER2 2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "SKBR3_HER2_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-HER2 3",
        name="SKBR3 HER2 3",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "SKBR3_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-PDL1 1",
        name="SKBR3 PD-L1 1",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "SKBR3_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-PDL1 2",
        name="SKBR3 PD-L1 2",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "SKBR3_PDL1_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-PDL1 3",
        name="SKBR3 PD-L1 3",
        molecule_conjugates=[("MBA", "PD-L1")],
        output_dir=output,
    ),
    "SKBR3_TROP2_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-TROP2 1",
        name="SKBR3 TROP2 1",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "SKBR3_TROP2_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-TROP2 2",
        name="SKBR3 TROP2 2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "SKBR3_TROP2_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-TROP2 3",
        name="SKBR3 TROP2 3",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "SKBR3_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-IgG 1",
        name="SKBR3 IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "SKBR3_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-IgG 2",
        name="SKBR3 IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "SKBR3_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3/SKBR3 - MBA-IgG 3",
        name="SKBR3 IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
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
        "mcf10a_epcam": [
            "MCF10A_EpCAM_1", "MCF10A_EpCAM_2", "MCF10A_EpCAM_3",
            "MCF10A_IgG_1", "MCF10A_IgG_2", "MCF10A_IgG_3",
        ],
        "mcf10a_her2": [
            "MCF10A_HER2_1", "MCF10A_HER2_2", "MCF10A_HER2_3",
            "MCF10A_IgG_1", "MCF10A_IgG_2", "MCF10A_IgG_3",
        ],
        "mcf10a_pdl1": [
            "MCF10A_PDL1_1", "MCF10A_PDL1_2", "MCF10A_PDL1_3",
            "MCF10A_IgG_1", "MCF10A_IgG_2", "MCF10A_IgG_3",
        ],
        "mcf10a_trop2": [
            "MCF10A_TROP2_1", "MCF10A_TROP2_2", "MCF10A_TROP2_3",
            "MCF10A_IgG_1", "MCF10A_IgG_2", "MCF10A_IgG_3",
        ],
        "mcf7_epcam": [
            "MCF7_EpCAM_1", "MCF7_EpCAM_2", "MCF7_EpCAM_3",
            "MCF7_IgG_1", "MCF7_IgG_2", "MCF7_IgG_3",
        ],
        "mcf7_her2": [
            "MCF7_HER2_1", "MCF7_HER2_2", "MCF7_HER2_3",
            "MCF7_IgG_1", "MCF7_IgG_2", "MCF7_IgG_3",
        ],
        "mcf7_pdl1": [
            "MCF7_PDL1_1", "MCF7_PDL1_2", "MCF7_PDL1_3",
            "MCF7_IgG_1", "MCF7_IgG_2", "MCF7_IgG_3",
        ],
        "mcf7_trop2": [
            "MCF7_TROP2_1", "MCF7_TROP2_2", "MCF7_TROP2_3",
            "MCF7_IgG_1", "MCF7_IgG_2", "MCF7_IgG_3",
        ],
        "skbr3_epcam": [
            "SKBR3_EpCAM_1", "SKBR3_EpCAM_2", "SKBR3_EpCAM_3",
            "SKBR3_IgG_1", "SKBR3_IgG_2", "SKBR3_IgG_3",
        ],
        "skbr3_her2": [
            "SKBR3_HER2_1", "SKBR3_HER2_2", "SKBR3_HER2_3",
            "SKBR3_IgG_1", "SKBR3_IgG_2", "SKBR3_IgG_3",
        ],
        "skbr3_pdl1": [
            "SKBR3_PDL1_1", "SKBR3_PDL1_2", "SKBR3_PDL1_3",
            "SKBR3_IgG_1", "SKBR3_IgG_2", "SKBR3_IgG_3",
        ],
        "skbr3_trop2": [
            "SKBR3_TROP2_1", "SKBR3_TROP2_2", "SKBR3_TROP2_3",
            "SKBR3_IgG_1", "SKBR3_IgG_2", "SKBR3_IgG_3",
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
