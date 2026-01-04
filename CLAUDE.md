# Raman Analysis Library - Context for Claude

## Project Overview

A Python library for analyzing Raman spectroscopy data with a focus on SERS (Surface Enhanced Raman Scattering) multiplex deconvolution. Used for cancer biomarker detection by analyzing multiplex antibody-conjugated nanoparticles.

**Domain context**: SERS uses Raman-active molecules (reporter molecules) attached to gold nanoparticles. Each molecule has a unique spectral fingerprint (characteristic peak). Multiplex samples contain multiple different molecules, and deconvolution determines the contribution of each.

## Core Design Principles

1. **Library over framework** - Provide composable functions, not a rigid workflow
2. **Python scripts over config files** - Experiment logic in Python, leveraging the library
3. **Automatic output management** - Create versioned output directories automatically
4. **Process all replicates** - Deconvolve each replicate individually for statistical analysis; use averaged spectra for visualization
5. **Minimal serialization** - Pure Python data structures (lists, dicts), no custom classes

## Raman Molecules

Three Raman-active molecules are currently supported:

```python
RAMAN_MOLECULES = {
    "MBA": {"peak": 1078, "color": "blue"},      # Mercaptobenzoic acid
    "DTNB": {"peak": 1335, "color": "green"},    # 5,5'-Dithiobis(2-nitrobenzoic acid)
    "TFMBA": {"peak": 1377, "color": "orange"}   # 4-(Trifluoromethyl)benzenethiol
}
```

Each has:
- **peak**: Characteristic Raman peak wavenumber (cm⁻¹)
- **color**: Standard color for plotting

## Analysis Workflow

### 1. Reference Samples (Pure single-molecule samples)
- Input: Directory with spectrum files, molecule tag (e.g., "MBA"), conjugate (e.g., "EpCAM")
- Process: Load → ALS baseline correction → Average → Plot
- Output: `reference_{molecule}_{conjugate}.png` showing raw and corrected spectra with expected peak

### 2. Experiment Samples (Multiplexed samples)
- Input: Directory with spectrum files, sample name, molecule-conjugate pairs list
- Process: Load → ALS baseline correction → Average → Store replicates → Plot
- Output: `sample_{sample_name}.png` showing averaged spectra with all expected peaks

### 3. Normalization & Deconvolution (Per Replicate)
- Input: All sample replicates + references, wavenumber range (typically 1000-1500 cm⁻¹)
- Process:
  - Normalize each replicate using L2 normalization
  - Deconvolve each replicate using NNLS to determine molecular contributions
  - Calculate mean ± std across all replicates
- Output:
  - `normalization_averaged_{sample}.png` - visualization using averaged spectrum
  - `deconvolution_averaged_{sample}.png` - 3-panel plot (multiplex vs fitted, individual contributions, residuals) in normalized scale
  - `deconvolution_averaged_original_scale_{sample}.png` - 3-panel plot in original intensity scale (before normalization)
  - `deconvolution_boxplots.png` - box plots showing distribution of contributions across all replicates for all samples

## Data Structures

All functions use simple Python dicts:

```python
# Single spectrum (after baseline correction)
spectrum = {
    "wavenumbers": [200.0, 201.0, ..., 2000.0],
    "intensities": [...],  # original
    "corrected": [...],    # baseline-corrected
    "baseline": [...]      # extracted baseline
}

# Averaged spectrum (for references)
averaged_reference = {
    "wavenumbers": [...],
    "raw_avg": [...],
    "raw_std": [...],
    "corrected_avg": [...],
    "corrected_std": [...],
    "baseline_avg": [...],
    "count": 150,  # number of spectra
    "molecule": "MBA",
    "conjugate": "EpCAM"
}

# Averaged spectrum (for samples) - includes replicates
averaged_sample = {
    "wavenumbers": [...],
    "raw_avg": [...],
    "raw_std": [...],
    "corrected_avg": [...],
    "corrected_std": [...],
    "baseline_avg": [...],
    "count": 150,  # number of spectra
    "molecule_conjugates": [("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
    "name": "Sample 1",
    "replicates": [  # List of individual corrected spectra
        {
            "wavenumbers": [...],
            "intensities": [...],
            "corrected": [...],
            "baseline": [...]
        },
        # ... 149 more replicates
    ]
}

# Deconvolution results (list of results, one per replicate)
deconv_results_for_sample = [
    {
        "contributions": {("MBA", "EpCAM"): 38.9, ("DTNB", "HER2"): 53.5, ("TFMBA", "TROP2"): 7.6},  # percentages
        "coefficients": {...},  # raw NNLS coefficients
        "reconstructed": [...],  # fitted spectrum (normalized)
        "residual": [...],
        "metrics": {"rmse": 3.09, "r_squared": 0.95},
        "individual_contributions": {("MBA", "EpCAM"): [...], ("DTNB", "HER2"): [...], ...},  # normalized
        "norm_factor": 123.45  # L2 norm used for sample normalization
    },
    # ... one result per replicate (150 total)
]
```

## Library Structure

```
raman_lib/
├── molecules.py       # RAMAN_MOLECULES config, get_peak(), get_color()
├── io.py              # load_spectra(), create_output_dir(), ensure_output_subdir()
├── baseline.py        # als_baseline(), apply_baseline_correction()
├── averaging.py       # calculate_average()
├── normalization.py   # normalize_l2(), normalize_spectra_l2()
├── deconvolution.py   # deconvolve_nnls()
├── plotting.py        # plot_reference(), plot_sample(), plot_normalization(), plot_deconvolution(), plot_deconvolution_original_scale(), plot_deconvolution_boxplots()
└── workflow.py        # build_reference_dict(), load_and_process_reference(), load_and_process_sample(), normalize_and_deconvolve_samples(), print_experiment_summary()
```

## Key API - Workflow Functions

These high-level functions encapsulate complete workflows:

```python
# Load and process references (pure single-molecule samples)
references = build_reference_dict([
    load_and_process_reference(
        directory="path/to/MBA_EpCAM",
        molecule="MBA",
        conjugate="EpCAM",
        output_dir=output
    ),
    load_and_process_reference(
        directory="path/to/DTNB_HER2",
        molecule="DTNB",
        conjugate="HER2",
        output_dir=output
    ),
    load_and_process_reference(
        directory="path/to/TFMBA_TROP2",
        molecule="TFMBA",
        conjugate="TROP2",
        output_dir=output
    ),
])

# Load and process a sample (multiplex)
samples = {
    "Sample_1": load_and_process_sample(
        directory="path/to/sample",
        name="Sample 1",
        molecule_conjugates=[("MBA", "EpCAM"), ("DTNB", "HER2"), ("TFMBA", "TROP2")],
        output_dir=output
    )
}

# Normalize and deconvolve all samples (processes each replicate individually)
deconv_results = normalize_and_deconvolve_samples(
    samples=samples,
    references=references,
    wavenumber_range=(1000, 1500),
    output_dir=output
)

# Print summary table
print_experiment_summary(
    output_dir=output,
    references=references,
    samples=samples,
    deconv_results=deconv_results
)
```

## Experiment Script Pattern

All experiment scripts follow this structure:

```python
#!/usr/bin/env python3
"""
Experiment: [Description]
Date: YYYY-MM-DD
"""

import os
from raman_lib import (
    create_output_dir,
    build_reference_dict,
    load_and_process_reference,
    load_and_process_sample,
    normalize_and_deconvolve_samples,
    print_experiment_summary
)

# Data directories
REFERENCE_DIR = os.path.expanduser("~/Documents/...")
SAMPLE_DIR = os.path.expanduser("~/Documents/...")
WAVENUMBER_RANGE = (1000, 1500)

# Create output (auto-versioned)
output = create_output_dir("experiment-name", base_dir="results")

# Load references (with explicit header in script)
print("="*60)
print("LOADING REFERENCES")
print("="*60)
references = {...}

# Load samples (with explicit header in script)
print("="*60)
print("LOADING SAMPLES")
print("="*60)
samples = {...}

# Normalize & deconvolve (with explicit header in script)
print("\n" + "="*60)
print("NORMALIZATION & DECONVOLUTION")
print("="*60)
deconv_results = normalize_and_deconvolve_samples(...)

# Summary
print_experiment_summary(...)
```

## Key Algorithms

- **ALS baseline**: Asymmetric Least Squares (Eilers & Boelens 2005). Iteratively fits smooth baseline by penalizing roughness and asymmetrically weighting points above/below curve.
  - lambda_param (1e7): controls smoothness
  - p (0.01): asymmetry (lower = more weight below baseline)
  - d (2): order of differences (second-order)

- **L2 normalization**: Divides spectrum by its Euclidean norm within specified wavenumber range. Makes spectra comparable regardless of intensity.

- **NNLS deconvolution**: Non-Negative Least Squares finds coefficients x ≥ 0 minimizing ||Ax - b||² where:
  - A = reference matrix (each column is a reference spectrum)
  - b = sample spectrum
  - x = contribution coefficients (converted to percentages)

## File Formats

- **Input spectra**: Tab-delimited `.txt` files with two columns (wavenumber, intensity)
- **Output plots**: PNG images at 300 DPI
- **Output directories**: Auto-versioned (name, name-2, name-3, etc.)
- **Output structure**: Flat directory with prefixed filenames:
  - `reference_{molecule}_{conjugate}.png` - Reference spectrum plots
  - `sample_{name}.png` - Sample spectrum plots (averaged)
  - `normalization_averaged_{name}.png` - Normalization comparison (averaged spectrum)
  - `deconvolution_averaged_{name}.png` - Deconvolution 3-panel plot (averaged spectrum, normalized scale)
  - `deconvolution_averaged_original_scale_{name}.png` - Deconvolution 3-panel plot (averaged spectrum, original scale)
  - `deconvolution_boxplots.png` - Box plots showing distribution across all replicates
  - `peak_intensity_histogram_{group}_{molecule}.png` - Peak intensity histograms per group/molecule

## Common Patterns

- Experiment files use date-prefixed naming: `YYYY_MM_DD_description.py`
- Output directories use date-prefixed naming: `YYYY-MM-DD-description`
- All directories created in `results/` (gitignored)
- References reused across experiments (no caching, point to raw data)
- Headers (like "LOADING REFERENCES") belong in experiment scripts, not library functions
- Workflow functions print detailed progress but not section summaries
- Individual replicates are processed for statistical analysis (mean ± std)
- Averaged spectra are used for visualization plots to verify algorithm performance
- Box plots show distribution of molecular contributions across all replicates
- References use `(molecule, conjugate)` tuples as dict keys, allowing multiple conjugates per molecule in one dict
- Samples specify their `molecule_conjugates` list, which is used to filter references for deconvolution
- Summary tables are grouped by conjugate type (e.g., all EpCAM/HER2/TROP2 samples in one table)
- References and samples are sorted alphabetically in output for easier reading
