#!/usr/bin/env python3
"""
Experiment: MDA-MB-231 and MDA-MB-468 Specificity - TFMBA-EpCAM and DTNB-PD-L1
Date: 2026-07-03
Description: Singleplex TFMBA-EpCAM and DTNB-PD-L1 peak intensity analysis of
             MDA-MB-231 and MDA-MB-468 cells with matched IgG controls.
             Histograms compare each cell line and nano tag against its matched
             IgG controls, producing four histograms total.
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
    "2026-07-03 MDA-MB-231 MDA-MB-468 - TFMBA-EpCAM DTNB-PDL1 IgG"
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
    # MDA-MB-231 TFMBA-EpCAM
    "MDA231_TFMBA_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 TFMA-EpCAM 1",
        name="MDA-MB-231 TFMBA EpCAM 1",
        molecule_conjugates=[("TFMBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_TFMBA_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 TFMA-EpCAM 2",
        name="MDA-MB-231 TFMBA EpCAM 2",
        molecule_conjugates=[("TFMBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_TFMBA_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 TFMBA-IgG 1",
        name="MDA-MB-231 TFMBA IgG 1",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_TFMBA_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 TFMBA-IgG 2",
        name="MDA-MB-231 TFMBA IgG 2",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    # MDA-MB-231 DTNB-PD-L1
    "MDA231_DTNB_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 DTNB-PDL1 1",
        name="MDA-MB-231 DTNB PD-L1 1",
        molecule_conjugates=[("DTNB", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_DTNB_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 DTNB-PDL1 2",
        name="MDA-MB-231 DTNB PD-L1 2",
        molecule_conjugates=[("DTNB", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_DTNB_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 DTNB-IgG 1",
        name="MDA-MB-231 DTNB IgG 1",
        molecule_conjugates=[("DTNB", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA231_DTNB_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-231 DTNB-IgG 2",
        name="MDA-MB-231 DTNB IgG 2",
        molecule_conjugates=[("DTNB", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    # MDA-MB-468 TFMBA-EpCAM
    "MDA468_TFMBA_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 TFMBA-EpCAM 1",
        name="MDA-MB-468 TFMBA EpCAM 1",
        molecule_conjugates=[("TFMBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_TFMBA_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 TFMBA-EpCAM 2",
        name="MDA-MB-468 TFMBA EpCAM 2",
        molecule_conjugates=[("TFMBA", "EpCAM")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_TFMBA_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 TFMBA-IgG 1",
        name="MDA-MB-468 TFMBA IgG 1",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_TFMBA_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 TFMBA-IgG 2",
        name="MDA-MB-468 TFMBA IgG 2",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    # MDA-MB-468 DTNB-PD-L1
    "MDA468_DTNB_PDL1_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 DTNB-PDL1 1",
        name="MDA-MB-468 DTNB PD-L1 1",
        molecule_conjugates=[("DTNB", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_DTNB_PDL1_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 DTNB-PDL1 2",
        name="MDA-MB-468 DTNB PD-L1 2",
        molecule_conjugates=[("DTNB", "PD-L1")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_DTNB_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 DTNB-IgG 1",
        name="MDA-MB-468 DTNB IgG 1",
        molecule_conjugates=[("DTNB", "IgG")],
        output_dir=output,
        raw_output_dir=RAW_OUTPUT_DIR,
    ),
    "MDA468_DTNB_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/MDA-MB-468 DTNB-IgG 2",
        name="MDA-MB-468 DTNB IgG 2",
        molecule_conjugates=[("DTNB", "IgG")],
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
        "mda231_tfmba_epcam": [
            "MDA231_TFMBA_EpCAM_1",
            "MDA231_TFMBA_EpCAM_2",
            "MDA231_TFMBA_IgG_1",
            "MDA231_TFMBA_IgG_2",
        ],
        "mda468_tfmba_epcam": [
            "MDA468_TFMBA_EpCAM_1",
            "MDA468_TFMBA_EpCAM_2",
            "MDA468_TFMBA_IgG_1",
            "MDA468_TFMBA_IgG_2",
        ],
    },
    molecules=["TFMBA"],
    x_max=325,
    output_dir=output,
    title_prefix="TFMBA-EpCAM",
)

plot_peak_histograms_from_samples(
    samples,
    groups={
        "mda231_dtnb_pdl1": [
            "MDA231_DTNB_PDL1_1",
            "MDA231_DTNB_PDL1_2",
            "MDA231_DTNB_IgG_1",
            "MDA231_DTNB_IgG_2",
        ],
        "mda468_dtnb_pdl1": [
            "MDA468_DTNB_PDL1_1",
            "MDA468_DTNB_PDL1_2",
            "MDA468_DTNB_IgG_1",
            "MDA468_DTNB_IgG_2",
        ],
    },
    molecules=["DTNB"],
    x_max=325,
    output_dir=output,
    title_prefix="DTNB-PD-L1",
)


# ============================================================
# Summary
# ============================================================

experiment_summary(
    experiment,
    samples=samples,
    output_dir=output,
)
