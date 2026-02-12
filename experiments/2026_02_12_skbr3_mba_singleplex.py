#!/usr/bin/env python3
"""
Experiment: SKBR3 MBA Singleplex — EpCAM, HER2, TROP2
Date: 2026-02-12
Description: Singleplex MBA peak intensity analysis across three antibody
             conjugates (EpCAM, HER2, TROP2) with BSA and IgG controls.
             All samples use MBA as the sole Raman reporter.
"""

import os
from raman_lib import (
    create_output_dir,
    load_and_process_sample,
    plot_peak_histograms_from_samples,
)

# Data directory
DATA_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/"
    "2026-02-12 SKBR3 MBA-EpCAM MBA-HER2 MBA-TROP2 Ag-coated AuNPs"
)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2026-02-12-skbr3-mba-singleplex", base_dir="results")
print(f"Output directory: {output}\n")


# ============================================================
# Step 1: Load SKBR3 cell samples
# ============================================================

print("=" * 60)
print("LOADING SKBR3 SAMPLES")
print("=" * 60)

skbr3_samples = {
    # EpCAM
    "MBA_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA EpCAM 1",
        name="MBA EpCAM 1",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA EpCAM 2",
        name="MBA EpCAM 2",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "MBA_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA EpCAM 3",
        name="MBA EpCAM 3",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    # HER2
    "MBA_HER2_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA HER2 1",
        name="MBA HER2 1",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MBA_HER2_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA HER2 2",
        name="MBA HER2 2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MBA_HER2_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA HER2 3",
        name="MBA HER2 3",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    # TROP2
    "MBA_TROP2_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA TROP2 1",
        name="MBA TROP2 1",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MBA_TROP2_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA TROP2 2",
        name="MBA TROP2 2",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "MBA_TROP2_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA TROP2 3",
        name="MBA TROP2 3",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    # BSA controls
    "MBA_BSA_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA BSA 1",
        name="MBA BSA 1",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_BSA_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA BSA 2",
        name="MBA BSA 2",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_BSA_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA BSA 3",
        name="MBA BSA 3",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    # IgG controls
    "MBA_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA IgG 1",
        name="MBA IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MBA_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA IgG 2",
        name="MBA IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MBA_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA IgG 3",
        name="MBA IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
}


# ============================================================
# Step 2: Load tags-only samples
# ============================================================

print("\n" + "=" * 60)
print("LOADING TAGS-ONLY SAMPLES")
print("=" * 60)

tag_samples = {
    "Tag_EpCAM": load_and_process_sample(
        f"{DATA_DIR}/MBA EpCAM tags",
        name="MBA EpCAM tags",
        molecule_conjugates=[("MBA", "EpCAM")],
        output_dir=output,
    ),
    "Tag_HER2": load_and_process_sample(
        f"{DATA_DIR}/MBA HER2 tags",
        name="MBA HER2 tags",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "Tag_TROP2": load_and_process_sample(
        f"{DATA_DIR}/MBA TROP2 tags",
        name="MBA TROP2 tags",
        molecule_conjugates=[("MBA", "TROP2")],
        output_dir=output,
    ),
    "Tag_BSA": load_and_process_sample(
        f"{DATA_DIR}/MBA BSA tags",
        name="MBA BSA tags",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "Tag_IgG": load_and_process_sample(
        f"{DATA_DIR}/MBA IgG tags",
        name="MBA IgG tags",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
}


# ============================================================
# Step 3: SKBR3 peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("SKBR3 PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_samples(
    skbr3_samples,
    groups={
        "epcam": [
            "MBA_EpCAM_1", "MBA_EpCAM_2", "MBA_EpCAM_3",
            "MBA_IgG_1", "MBA_IgG_2", "MBA_IgG_3",
            "MBA_BSA_1", "MBA_BSA_2", "MBA_BSA_3",
        ],
        "her2": [
            "MBA_HER2_1", "MBA_HER2_2", "MBA_HER2_3",
            "MBA_IgG_1", "MBA_IgG_2", "MBA_IgG_3",
            "MBA_BSA_1", "MBA_BSA_2", "MBA_BSA_3",
        ],
        "trop2": [
            "MBA_TROP2_1", "MBA_TROP2_2", "MBA_TROP2_3",
            "MBA_IgG_1", "MBA_IgG_2", "MBA_IgG_3",
            "MBA_BSA_1", "MBA_BSA_2", "MBA_BSA_3",
        ],
    },
    molecules=["MBA"],
    output_dir=output,
)


# ============================================================
# Step 4: Tags-only peak intensity histogram
# ============================================================

print("\n" + "=" * 60)
print("TAGS-ONLY PEAK INTENSITY HISTOGRAM")
print("=" * 60)

plot_peak_histograms_from_samples(
    tag_samples,
    groups={
        "tags": ["Tag_EpCAM", "Tag_HER2", "Tag_TROP2", "Tag_IgG", "Tag_BSA"],
    },
    molecules=["MBA"],
    output_dir=output,
    bin_size=250
)


# ============================================================
# Summary
# ============================================================

print("\n" + "=" * 60)
print("EXPERIMENT COMPLETE")
print("=" * 60)

print(f"\nOutput directory: {output}")
print(f"\nSKBR3 samples processed:")
for sample_key, data in sorted(skbr3_samples.items(), key=lambda x: x[1]["name"]):
    print(f"  {data['name']}: {data['count']} spectra")
print(f"\nTag samples processed:")
for sample_key, data in sorted(tag_samples.items(), key=lambda x: x[1]["name"]):
    print(f"  {data['name']}: {data['count']} spectra")
