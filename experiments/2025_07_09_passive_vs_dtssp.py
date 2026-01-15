#!/usr/bin/env python3
"""
Experiment: Passive vs DTSSP Antibody Conjugation Comparison
Date: 2025-07-09
Description: Comparing passive vs DTSSP antibody conjugation methods using MBA-EpCAM on MCF7 cells.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
)

# Data directories
DTSSP_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-07-09 passive vs DTSSP Ab conjugation, MBA EpCAM, MCF7/DTSSP conjugation"
)
PASSIVE_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/2025-07-09 passive vs DTSSP Ab conjugation, MBA EpCAM, MCF7/Passive Conjugation"
)

# Create output directory
output = create_output_dir("2025-07-09-passive-vs-dtssp", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Load and process samples
# ============================================================

print("=" * 60)
print("LOADING SAMPLES")
print("=" * 60)

samples = {
    # DTSSP conjugation
    "DTSSP_EpCAM_2": load_and_process_sample(
        f"{DTSSP_DIR}/MBA EpCAM DTSSP 2",
        name="DTSSP EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "DTSSP_IgG_1": load_and_process_sample(
        f"{DTSSP_DIR}/MBA IgG DTSSP 1",
        name="DTSSP IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "DTSSP_IgG_2": load_and_process_sample(
        f"{DTSSP_DIR}/MBA IgG DTSSP 2",
        name="DTSSP IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "DTSSP_BSA_1": load_and_process_sample(
        f"{DTSSP_DIR}/MBA DTSSP BSA 1 mistaken EpCAM maybe passive",
        name="DTSSP BSA 1 (mistaken)",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "DTSSP_BSA_2": load_and_process_sample(
        f"{DTSSP_DIR}/MBA DTSSP BSA 2 mistaken EpCAM maybe passive",
        name="DTSSP BSA 2 (mistaken)",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    # Passive conjugation
    "Passive_EpCAM_2": load_and_process_sample(
        f"{PASSIVE_DIR}/MBA EpCAM 2 passive",
        name="Passive EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Passive_IgG_2": load_and_process_sample(
        f"{PASSIVE_DIR}/MBA IgG 2 passive",
        name="Passive IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
}

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
