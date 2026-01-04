# Raman Analysis Library

A Python library for analyzing Raman spectroscopy data with a focus on SERS (Surface Enhanced Raman Scattering) multiplex deconvolution.

## Installation

```bash
uv sync
```

## Usage

Create an experiment script that processes reference spectra and multiplex samples:

```python
#!/usr/bin/env python3
"""
Experiment: SKBR3 Spiked PBMCs Multiplex Analysis
Date: 2025-11-06
"""

import os
from raman_lib import (
    create_output_dir,
    build_reference_dict,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples,
    plot_peak_histograms_from_deconv,
    print_experiment_summary,
)

# Data directories
REFERENCE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-07-17 SKBR3 passive")
SAMPLE_DIR = os.path.expanduser("~/Documents/Spectroscopy Results/2025-11-06 SKBR3 spiked PBMCs multiplex")

# Wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Create output directory (auto-versioned) in results/
output = create_output_dir("2025-11-06-skbr3-spiked-pbmc", base_dir="results")
print(f"Output directory: {output}\n")

# ============================================================
# Step 1: Load and process references (from 2025-07-17)
# ============================================================

print("="*60)
print("LOADING REFERENCES")
print("="*60)

references = build_reference_dict([
    load_and_process_reference(
        f"{REFERENCE_DIR}/MBA EpCAM",
        molecule="MBA",
        conjugate="EpCAM",
        output_dir=output
    ),
    load_and_process_reference(
        f"{REFERENCE_DIR}/DTNB HER2",
        molecule="DTNB",
        conjugate="HER2",
        output_dir=output
    ),
    load_and_process_reference(
        f"{REFERENCE_DIR}/TFMBA TROP2",
        molecule="TFMBA",
        conjugate="TROP2",
        output_dir=output
    ),
])


# ============================================================
# Step 2: Load and process multiplex samples (from 2025-11-06)
# ============================================================

print("="*60)
print("LOADING SAMPLES")
print("="*60)

samples = {
    "Multiplex_Ab_1": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 1",
        name="Multiplex Ab 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output
    ),
    "Multiplex_Ab_2": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 2",
        name="Multiplex Ab 2",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output
    ),
    "Multiplex_Ab_3": load_and_process_sample(
        f"{SAMPLE_DIR}/SKBR3 spiked PBMC_10E5_Multiplex Ab 3",
        name="Multiplex Ab 3",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output
    ),
}


# ============================================================
# Step 3: Normalize and deconvolve all replicates
# ============================================================

print("\n" + "="*60)
print("NORMALIZATION & DECONVOLUTION")
print("="*60)

deconv_results = normalize_and_deconvolve_samples(
    samples=samples,
    references=references,
    wavenumber_range=WAVENUMBER_RANGE,
    output_dir=output
)


# ============================================================
# Peak intensity histograms
# ============================================================

print("\n" + "="*60)
print("PEAK INTENSITY HISTOGRAMS")
print("="*60)

plot_peak_histograms_from_deconv(
    deconv_results,
    groups={"multiplex": ["Multiplex_Ab_1", "Multiplex_Ab_2", "Multiplex_Ab_3"]},
    output_dir=output,
)


# ============================================================
# Summary
# ============================================================

print("\n" + "="*60)
print("EXPERIMENT COMPLETE")
print("="*60)

print_experiment_summary(
    output_dir=output,
    references=references,
    samples=samples,
    deconv_results=deconv_results
)
```

Run your experiment:

```bash
uv run python experiments/2025_11_06_skbr3_spiked_pbmc.py
```

## Output

Results are saved in auto-versioned directories under `results/` with a flat structure:

```
results/2025-11-06-skbr3-spiked-pbmc/
├── reference_MBA_EpCAM.png
├── reference_DTNB_HER2.png
├── reference_TFMBA_TROP2.png
├── sample_Multiplex_Ab_1.png
├── sample_Multiplex_Ab_2.png
├── sample_Multiplex_Ab_3.png
├── normalization_averaged_Multiplex_Ab_1.png
├── normalization_averaged_Multiplex_Ab_2.png
├── normalization_averaged_Multiplex_Ab_3.png
├── deconvolution_averaged_Multiplex_Ab_1.png
├── deconvolution_averaged_Multiplex_Ab_2.png
├── deconvolution_averaged_Multiplex_Ab_3.png
├── deconvolution_averaged_original_scale_Multiplex_Ab_1.png
├── deconvolution_averaged_original_scale_Multiplex_Ab_2.png
├── deconvolution_averaged_original_scale_Multiplex_Ab_3.png
├── deconvolution_boxplots.png
├── peak_intensity_histogram_multiplex_MBA.png
├── peak_intensity_histogram_multiplex_DTNB.png
└── peak_intensity_histogram_multiplex_TFMBA.png
```

The workflow performs:
1. **Baseline correction** (ALS algorithm) on all spectra
2. **Averaging** of replicate spectra (for visualization)
3. **L2 normalization** of each individual replicate within specified wavenumber range
4. **NNLS deconvolution** on each replicate to determine molecular contributions
5. **Statistical analysis** calculating mean ± standard deviation across all replicates
6. **Box plots** showing the distribution of contributions for quality assessment
7. **Peak intensity histograms** showing intensity distributions at characteristic peaks
