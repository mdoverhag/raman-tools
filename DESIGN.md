# Raman Analysis Library - Design Document

## Overview

A Python library for analyzing Raman spectroscopy data with a focus on SERS (Surface Enhanced Raman Scattering) multiplex deconvolution. The goal is to provide reusable algorithms that can be composed in experiment scripts without the rigidity of a GUI application or the manual complexity of separate command-line tools.

## Key Design Principles

1. **Library over framework** - Provide composable functions, not a rigid workflow
2. **Python scripts over config files** - Experiment logic in Python, leveraging the library
3. **Automatic output management** - Create versioned output directories automatically
4. **Work with averages first** - Process averaged spectra, not individual replicates
5. **Minimal serialization** - Pure Python data structures (lists, dicts), no custom classes

## Molecule Configuration

```python
# Built into the library
RAMAN_MOLECULES = {
    "MBA": {"peak": 1078, "color": "blue"},
    "DTNB": {"peak": 1335, "color": "green"},
    "TFMBA": {"peak": 1377, "color": "orange"}
}
```

## Workflow

### Step 1: Reference Samples (Pure single-molecule samples)

**Input**:
- Directory path to reference spectrum files
- Molecule tag (e.g., "MBA", "DTNB", "TFMBA")

**Process**:
1. Load all spectra in directory
2. Apply ALS baseline correction
3. Calculate average of raw spectra
4. Calculate average of baseline-corrected spectra

**Output**: `references/{molecule_name}.png` - Plot showing:
- Average raw spectrum
- Average corrected spectrum
- Expected peak location highlighted (using molecule tag)

### Step 2: Experiment Samples (Multiplexed samples)

**Input**:
- Directory path to sample spectrum files
- Sample name (e.g., "Multi Ab 1", "Multi Ab 2", "Multi Ab 3")
- Molecule tags list (typically ["MBA", "DTNB", "TFMBA"] for multiplex)

**Process**:
1. Load all spectra in directory
2. Apply ALS baseline correction
3. Calculate average of raw spectra
4. Calculate average of baseline-corrected spectra

**Output**: `samples/{sample_name}.png` - Plot showing:
- Average raw spectrum
- Average corrected spectrum
- Expected peaks for all molecules in tags highlighted

### Step 3: Normalization

**Input**:
- Averaged experiment sample spectrum (with molecule tags)
- Averaged reference spectra (from Step 1)
- Normalization method ("l2" or "minmax")
- Wavenumber range for normalization (script-defined, e.g., (1000, 1500))

**Process**:
1. Apply normalization method (L2 or MinMax) to sample and all references
2. Focus on specified wavenumber range

**Output**: `deconvolution/normalization_{sample_name}.png` - Plot showing:
- Normalized sample spectrum
- Normalized reference spectra overlaid (only those in sample's molecule tags)
- Normalization range highlighted
- Expected peaks marked

### Step 4: Deconvolution (NNLS)

**Input**:
- Normalized sample spectrum
- Normalized reference spectra
- Wavenumber range for analysis

**Process**:
1. Perform Non-Negative Least Squares fitting
2. Calculate contribution percentages for each Raman molecule
3. Calculate reconstruction quality metrics (RMSE, R²)

**Output**: `deconvolution/deconvolution_{sample_name}.png` - Three-panel plot:
- **Top panel**: Multiplex spectrum (black) vs Fitted reconstruction (red dashed)
- **Middle panel**: Individual SERS tag contributions with percentages in legend
  - MBA in blue
  - DTNB in green
  - TFMBA in orange
- **Bottom panel**: Fitting residuals with RMS value

## Library Structure

```
raman_lib/
├── __init__.py
├── molecules.py        # RAMAN_MOLECULES configuration
├── io.py              # load_spectra(), save_plot()
├── baseline.py        # als_baseline() - ALS only
├── averaging.py       # calculate_average()
├── normalization.py   # normalize_l2(), normalize_minmax()
├── deconvolution.py   # nnls_deconvolve()
└── plotting.py        # plot_reference(), plot_sample(), plot_normalization(), plot_deconvolution()
```

## Experiment Script Template

```python
#!/usr/bin/env python3
"""
Experiment: [Description]
Date: YYYY-MM-DD
"""

from raman_lib import *

# Create output directory (auto-versioned: experiment-name, experiment-name-2, etc.)
output = create_output_dir("experiment-name")

# Define wavenumber range for normalization and deconvolution
WAVENUMBER_RANGE = (1000, 1500)

# Step 1: Load and process references (tag with molecule)
references = {
    "MBA": load_and_process_reference(
        "data/2025-01-15/MBA_EpCAM",
        molecule="MBA",
        output_dir=output
    ),
    "DTNB": load_and_process_reference(
        "data/2025-01-15/DTNB_HER2",
        molecule="DTNB",
        output_dir=output
    ),
    "TFMBA": load_and_process_reference(
        "data/2025-01-15/TFMBA_TROP2",
        molecule="TFMBA",
        output_dir=output
    ),
}

# Step 2: Load and process samples (tag with molecule list)
samples = {
    "Multi_Ab_1": load_and_process_sample(
        "data/2025-01-15/Multi_Ab_1",
        name="Multi Ab 1",
        molecules=["MBA", "DTNB", "TFMBA"],  # all three in multiplex
        output_dir=output
    ),
    "Multi_Ab_2": load_and_process_sample(
        "data/2025-01-15/Multi_Ab_2",
        name="Multi Ab 2",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
    "Multi_Ab_3": load_and_process_sample(
        "data/2025-01-15/Multi_Ab_3",
        name="Multi Ab 3",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_dir=output
    ),
}

# Step 3 & 4: Normalize and deconvolute each sample
for name, sample in samples.items():
    normalized = normalize_spectra(
        sample=sample,
        references=references,
        method="l2",  # or "minmax"
        wavenumber_range=WAVENUMBER_RANGE,
        output_dir=output
    )

    result = deconvolve_nnls(
        sample=normalized["sample"],
        references=normalized["references"],
        wavenumber_range=WAVENUMBER_RANGE,
        output_dir=output
    )

    # Plots saved automatically to output directory
```

## Output Directory Structure

```
experiment-name/
├── references/
│   ├── MBA.png
│   ├── DTNB.png
│   └── TFMBA.png
├── samples/
│   ├── Multi_Ab_1.png
│   ├── Multi_Ab_2.png
│   └── Multi_Ab_3.png
└── deconvolution/
    ├── normalization_Multi_Ab_1.png
    ├── deconvolution_Multi_Ab_1.png
    ├── normalization_Multi_Ab_2.png
    ├── deconvolution_Multi_Ab_2.png
    ├── normalization_Multi_Ab_3.png
    └── deconvolution_Multi_Ab_3.png
```

## Data Structures

All functions use simple Python types:

```python
# Spectrum data
spectrum = {
    "wavenumbers": [200.0, 201.0, ..., 2000.0],  # 1801 points
    "intensities": [123.4, 125.6, ..., 98.7],
}

# Averaged spectrum with statistics
averaged_spectrum = {
    "wavenumbers": [200.0, 201.0, ..., 2000.0],
    "raw_avg": [...],
    "raw_std": [...],
    "corrected_avg": [...],
    "corrected_std": [...],
    "baseline_avg": [...],
    "count": 150,  # number of spectra averaged
    "molecule": "MBA",  # for references
    # OR
    "molecules": ["MBA", "DTNB", "TFMBA"],  # for samples
}

# Deconvolution result
deconv_result = {
    "contributions": {"MBA": 38.9, "DTNB": 53.5, "TFMBA": 7.6},
    "reconstructed": [...],
    "residual": [...],
    "metrics": {"rmse": 3.09, "r_squared": 0.95},
}
```

## Key Library Functions (API Sketch)

```python
# High-level workflow functions
def load_and_process_reference(directory: str, molecule: str, output_dir: str) -> dict:
    """Load reference spectra, apply ALS baseline, average, and plot."""

def load_and_process_sample(directory: str, name: str, molecules: list[str], output_dir: str) -> dict:
    """Load sample spectra, apply ALS baseline, average, and plot."""

def normalize_spectra(sample: dict, references: dict, method: str, wavenumber_range: tuple, output_dir: str) -> dict:
    """Normalize sample and references, create plot."""

def deconvolve_nnls(sample: dict, references: dict, wavenumber_range: tuple, output_dir: str) -> dict:
    """Perform NNLS deconvolution, create plot."""

# Low-level utility functions
def load_spectra(directory: str) -> list[dict]:
    """Load all .txt files from directory."""

def als_baseline(spectrum: list[float], lambda_param: float, p: float, d: int) -> tuple[list, list]:
    """Apply ALS baseline correction."""

def calculate_average(spectra: list[dict]) -> dict:
    """Calculate mean and std of multiple spectra."""

def plot_deconvolution(sample, references, result, filepath: str):
    """Create 3-panel deconvolution plot."""
```

## Design Decisions (Resolved)

1. ✅ **Baseline correction method**: ALS only (no SNIP)
2. ✅ **Normalization methods**: L2 and MinMax (skip Area)
3. ✅ **Wavenumber range**: Configurable per experiment script
4. ✅ **File naming**: Simple descriptive names, no timestamps
5. ✅ **Reference reuse**: No caching, point to raw data in each script
6. ✅ **Molecule tagging**:
   - References: single molecule string
   - Samples: list of molecules present

## Migration from Previous Projects

**From spectometry-graphs**:
- ✅ Keep: ALS baseline, NNLS deconvolution algorithms
- ✅ Keep: Three-panel deconvolution plot style
- ❌ Drop: Manual multi-script workflow
- ❌ Drop: CSV intermediate files (work in memory)

**From raman-tools**:
- ✅ Keep: Python algorithm implementations (baseline, normalization, NNLS)
- ✅ Keep: Clean separation of concerns
- ❌ Drop: Rust/Tauri GUI
- ❌ Drop: JSON serialization overhead
- ❌ Drop: Complex state management

## Implementation Plan

1. **Phase 1**: Create library structure and basic I/O functions
2. **Phase 2**: Implement baseline correction and averaging
3. **Phase 3**: Implement normalization methods
4. **Phase 4**: Implement NNLS deconvolution
5. **Phase 5**: Create plotting functions
6. **Phase 6**: Build high-level workflow functions
7. **Phase 7**: Create experiment template and test with real data
