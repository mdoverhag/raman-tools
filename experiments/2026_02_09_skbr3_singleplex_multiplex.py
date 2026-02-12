#!/usr/bin/env python3
"""
Experiment: SKBR3 Singleplex and Multiplex — TFMBA-EpCAM, MBA-HER2, MMTAA-TROP2
Date: 2026-02-09
Description: Singleplex peak intensity analysis and multiplex deconvolution
             of SKBR3 cells with three Raman reporters. First experiment
             using MMTAA (peak 1287 cm⁻¹) replacing DTNB.
"""

import os
from raman_lib import (
    create_output_dir,
    build_reference_dict,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples,
    plot_peak_histograms_from_samples,
    plot_peak_histograms_from_deconv,
    print_experiment_summary,
)

# Data directories
DATA_DIR = os.path.expanduser(
    "~/Documents/Spectroscopy Results/"
    "2026-02-09 SKBR3 TFMBA-Ep MBA-H2 MMTAA T2 Singleplex and Multiplex"
)
REFERENCE_DIR = f"{DATA_DIR}/Ag-coated AuNP- Ab tags only"
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2026-02-09-skbr3-singleplex-multiplex", base_dir="results")
print(f"Output directory: {output}\n")


# ============================================================
# Step 1: Singleplex analysis — peak intensity histograms
# ============================================================

print("=" * 60)
print("LOADING SINGLEPLEX SAMPLES")
print("=" * 60)

singleplex_samples = {
    # MBA singleplex (HER2-conjugated)
    "MBA_HER2_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA H1",
        name="MBA HER2 1",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MBA_HER2_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA H2",
        name="MBA HER2 2",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MBA_HER2_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA H3",
        name="MBA HER2 3",
        molecule_conjugates=[("MBA", "HER2")],
        output_dir=output,
    ),
    "MBA_BSA_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA B1",
        name="MBA BSA 1",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_BSA_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA B2",
        name="MBA BSA 2",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_BSA_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA B3",
        name="MBA BSA 3",
        molecule_conjugates=[("MBA", "BSA")],
        output_dir=output,
    ),
    "MBA_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA G1",
        name="MBA IgG 1",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MBA_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA G2",
        name="MBA IgG 2",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    "MBA_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MBA G3",
        name="MBA IgG 3",
        molecule_conjugates=[("MBA", "IgG")],
        output_dir=output,
    ),
    # TFMBA singleplex (EpCAM-conjugated)
    "TFMBA_EpCAM_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA E1",
        name="TFMBA EpCAM 1",
        molecule_conjugates=[("TFMBA", "EpCAM")],
        output_dir=output,
    ),
    "TFMBA_EpCAM_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA E2",
        name="TFMBA EpCAM 2",
        molecule_conjugates=[("TFMBA", "EpCAM")],
        output_dir=output,
    ),
    "TFMBA_EpCAM_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA E3",
        name="TFMBA EpCAM 3",
        molecule_conjugates=[("TFMBA", "EpCAM")],
        output_dir=output,
    ),
    "TFMBA_BSA_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA B1",
        name="TFMBA BSA 1",
        molecule_conjugates=[("TFMBA", "BSA")],
        output_dir=output,
    ),
    "TFMBA_BSA_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA B2",
        name="TFMBA BSA 2",
        molecule_conjugates=[("TFMBA", "BSA")],
        output_dir=output,
    ),
    "TFMBA_BSA_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA B3",
        name="TFMBA BSA 3",
        molecule_conjugates=[("TFMBA", "BSA")],
        output_dir=output,
    ),
    "TFMBA_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA G1",
        name="TFMBA IgG 1",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
    ),
    "TFMBA_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA G2",
        name="TFMBA IgG 2",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
    ),
    "TFMBA_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 TFMBA G3",
        name="TFMBA IgG 3",
        molecule_conjugates=[("TFMBA", "IgG")],
        output_dir=output,
    ),
    # MMTAA singleplex (TROP2-conjugated)
    "MMTAA_TROP2_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA T1",
        name="MMTAA TROP2 1",
        molecule_conjugates=[("MMTAA", "TROP2")],
        output_dir=output,
    ),
    "MMTAA_TROP2_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA T2",
        name="MMTAA TROP2 2",
        molecule_conjugates=[("MMTAA", "TROP2")],
        output_dir=output,
    ),
    "MMTAA_TROP2_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA T3",
        name="MMTAA TROP2 3",
        molecule_conjugates=[("MMTAA", "TROP2")],
        output_dir=output,
    ),
    "MMTAA_BSA_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA B1",
        name="MMTAA BSA 1",
        molecule_conjugates=[("MMTAA", "BSA")],
        output_dir=output,
    ),
    "MMTAA_BSA_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA B2",
        name="MMTAA BSA 2",
        molecule_conjugates=[("MMTAA", "BSA")],
        output_dir=output,
    ),
    "MMTAA_BSA_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA B3",
        name="MMTAA BSA 3",
        molecule_conjugates=[("MMTAA", "BSA")],
        output_dir=output,
    ),
    "MMTAA_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA G1",
        name="MMTAA IgG 1",
        molecule_conjugates=[("MMTAA", "IgG")],
        output_dir=output,
    ),
    "MMTAA_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA G2",
        name="MMTAA IgG 2",
        molecule_conjugates=[("MMTAA", "IgG")],
        output_dir=output,
    ),
    "MMTAA_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 MMTAA G3",
        name="MMTAA IgG 3",
        molecule_conjugates=[("MMTAA", "IgG")],
        output_dir=output,
    ),
}

print("\n" + "=" * 60)
print("SINGLEPLEX PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_samples(
    singleplex_samples,
    groups={
        "singleplex": [
            "MBA_HER2_1", "MBA_HER2_2", "MBA_HER2_3",
            "MBA_BSA_1", "MBA_BSA_2", "MBA_BSA_3",
            "MBA_IgG_1", "MBA_IgG_2", "MBA_IgG_3",
        ],
    },
    output_dir=output,
    molecules=["MBA"],
)

plot_peak_histograms_from_samples(
    singleplex_samples,
    groups={
        "singleplex": [
            "TFMBA_EpCAM_1", "TFMBA_EpCAM_2", "TFMBA_EpCAM_3",
            "TFMBA_BSA_1", "TFMBA_BSA_2", "TFMBA_BSA_3",
            "TFMBA_IgG_1", "TFMBA_IgG_2", "TFMBA_IgG_3",
        ],
    },
    output_dir=output,
    molecules=["TFMBA"],
)

plot_peak_histograms_from_samples(
    singleplex_samples,
    groups={
        "singleplex": [
            "MMTAA_TROP2_1", "MMTAA_TROP2_2", "MMTAA_TROP2_3",
            "MMTAA_BSA_1", "MMTAA_BSA_2", "MMTAA_BSA_3",
            "MMTAA_IgG_1", "MMTAA_IgG_2", "MMTAA_IgG_3",
        ],
    },
    output_dir=output,
    molecules=["MMTAA"],
)


# ============================================================
# Step 2: Load and process references (from this experiment)
# ============================================================

print("\n" + "=" * 60)
print("LOADING REFERENCES")
print("=" * 60)

references = build_reference_dict(
    [
        # Antibody-conjugated references
        load_and_process_reference(
            f"{REFERENCE_DIR}/MBA HER2 tag",
            molecule="MBA",
            conjugate="HER2",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/TFMBA EpCAM tag",
            molecule="TFMBA",
            conjugate="EpCAM",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/MMTAA TROP2 tag",
            molecule="MMTAA",
            conjugate="TROP2",
            output_dir=output,
        ),
        # BSA control references
        load_and_process_reference(
            f"{REFERENCE_DIR}/MBA BSA tag",
            molecule="MBA",
            conjugate="BSA",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/TFMBA BSA tag",
            molecule="TFMBA",
            conjugate="BSA",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/MMTAA BSA tag",
            molecule="MMTAA",
            conjugate="BSA",
            output_dir=output,
        ),
        # IgG control references
        load_and_process_reference(
            f"{REFERENCE_DIR}/MBA IgG tag",
            molecule="MBA",
            conjugate="IgG",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/TFMBA IgG tag",
            molecule="TFMBA",
            conjugate="IgG",
            output_dir=output,
        ),
        load_and_process_reference(
            f"{REFERENCE_DIR}/MMTAA IgG tag",
            molecule="MMTAA",
            conjugate="IgG",
            output_dir=output,
        ),
    ]
)


# ============================================================
# Step 3: Load and process multiplex samples
# ============================================================

print("\n" + "=" * 60)
print("LOADING MULTIPLEX SAMPLES")
print("=" * 60)

multiplex_samples = {
    "Multi_Ab_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi Ab1",
        name="Multi Ab 1",
        molecule_conjugates=[("TFMBA", "EpCAM"), ("MBA", "HER2"), ("MMTAA", "TROP2")],
        output_dir=output,
    ),
    "Multi_Ab_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi Ab2",
        name="Multi Ab 2",
        molecule_conjugates=[("TFMBA", "EpCAM"), ("MBA", "HER2"), ("MMTAA", "TROP2")],
        output_dir=output,
    ),
    "Multi_Ab_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi Ab3",
        name="Multi Ab 3",
        molecule_conjugates=[("TFMBA", "EpCAM"), ("MBA", "HER2"), ("MMTAA", "TROP2")],
        output_dir=output,
    ),
    "Multi_BSA_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi B1",
        name="Multi BSA 1",
        molecule_conjugates=[("TFMBA", "BSA"), ("MBA", "BSA"), ("MMTAA", "BSA")],
        output_dir=output,
    ),
    "Multi_BSA_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi B2",
        name="Multi BSA 2",
        molecule_conjugates=[("TFMBA", "BSA"), ("MBA", "BSA"), ("MMTAA", "BSA")],
        output_dir=output,
    ),
    "Multi_BSA_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi B3",
        name="Multi BSA 3",
        molecule_conjugates=[("TFMBA", "BSA"), ("MBA", "BSA"), ("MMTAA", "BSA")],
        output_dir=output,
    ),
    "Multi_IgG_1": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi G1",
        name="Multi IgG 1",
        molecule_conjugates=[("TFMBA", "IgG"), ("MBA", "IgG"), ("MMTAA", "IgG")],
        output_dir=output,
    ),
    "Multi_IgG_2": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi G2",
        name="Multi IgG 2",
        molecule_conjugates=[("TFMBA", "IgG"), ("MBA", "IgG"), ("MMTAA", "IgG")],
        output_dir=output,
    ),
    "Multi_IgG_3": load_and_process_sample(
        f"{DATA_DIR}/SKBR3 Multi G3",
        name="Multi IgG 3",
        molecule_conjugates=[("TFMBA", "IgG"), ("MBA", "IgG"), ("MMTAA", "IgG")],
        output_dir=output,
    ),
}


# ============================================================
# Step 4: Normalize and deconvolve multiplex samples
# ============================================================

print("\n" + "=" * 60)
print("NORMALIZATION & DECONVOLUTION")
print("=" * 60)

deconv_results = normalize_and_deconvolve_samples(
    samples=multiplex_samples,
    references=references,
    wavenumber_range=WAVENUMBER_RANGE,
    output_dir=output,
)


# ============================================================
# Step 5: Multiplex peak intensity histograms
# ============================================================

print("\n" + "=" * 60)
print("MULTIPLEX PEAK INTENSITY HISTOGRAMS")
print("=" * 60)

plot_peak_histograms_from_deconv(
    deconv_results,
    groups={
        "multiplex": [
            "Multi_Ab_1", "Multi_Ab_2", "Multi_Ab_3",
            "Multi_BSA_1", "Multi_BSA_2", "Multi_BSA_3",
            "Multi_IgG_1", "Multi_IgG_2", "Multi_IgG_3",
        ],
    },
    output_dir=output,
    bin_size=10
)


# ============================================================
# Summary
# ============================================================

print("\n" + "=" * 60)
print("EXPERIMENT COMPLETE")
print("=" * 60)

print_experiment_summary(
    output_dir=output,
    references=references,
    samples=multiplex_samples,
    deconv_results=deconv_results,
)
