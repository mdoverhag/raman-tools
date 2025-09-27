# Changelog

All notable changes to Raman Tools will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2025-09-28

Adds an end-to-end NNLS deconvolution workflow (with L2 normalization) spanning new Python scripts, Rust storage/commands, and Svelte UI for visualization; also alphabetically sorts samples.

- **Backend**
  - **Deconvolution pipeline**: Add `DeconvolutionRun`/`DeconvolutionStorage` and Tauri commands: `create_deconvolution_run`, `get_deconvolution_run`, `calculate_normalization`, `perform_deconvolution` in `src-tauri/src/lib.rs` and `deconvolution.rs`.
  - **Python bridge**: Implement `normalize_spectra(...)` and `perform_nnls_deconvolution(...)` with request/response structs in `python_bridge.rs`.
  - **Python scripts**: Add `src-tauri/python/normalize_spectra.py` (L2 normalization in `1000–1500 cm⁻¹`) and `src-tauri/python/deconvolute_nnls.py` (NNLS using `scipy.optimize.nnls`).
  - **Python setup**: Expose `get_runtime_dir` and sync new scripts in `python_setup.rs`.

- **Frontend**
  - **UI components**: Add `DeconvolutionView.svelte` and `NormalizationChart.svelte` with normalized spectra, component contributions, reconstruction, and residual plots; validates required singleplex references; shows L2 norm table.
  - **Navigation**: Add "Deconvolution" tab in `+page.svelte`; extend store with `selectDeconvolution()`.
  - **Sorting**: Alphabetically sort samples in `SampleSidebar.svelte`.

## [0.5.0] - 2025-09-22

### Added

- Automatic average spectrum calculation after importing files
- Average spectrum displayed at top of spectra list with special styling
- Auto-selection of average spectrum for cleaner initial view
- Separate storage of raw and baseline-corrected averages

## [0.4.0] - 2025-09-14

### Added

- New `molecules.rs` module with `MoleculePair` struct for Raman-Target molecule pairs
- `MoleculePairEditor.svelte` component for editing molecule pairs in UI

### Changed

- Raman and Target molecules now stored as paired relationships instead of separate arrays
- Updated `Sample` struct to use `Vec<MoleculePair>` instead of separate `raman_molecules` and `target_molecules` fields
- Moved molecule pair editing UI from `+page.svelte` to dedicated `MoleculePairEditor.svelte` component
- Simplified `SampleSidebar.svelte` by removing inline molecule editing

### Removed

- `molecule-colors.ts` utility file (color logic moved to component)
- Separate `raman_molecules` and `target_molecules` fields from Sample struct

## [0.3.0] - 2025-09-10

### Added

- Batch processing for baseline correction - single Python subprocess handles all spectra
- New `batch_processor.py` script that streams JSON results as each spectrum completes
- `batch_baseline.rs` module providing iterator-based interface for batch processing
- `spectrum_importer.rs` module to handle file imports and emit progress events
- `samples.rs` module with CRUD operations for sample management
- Three-stage progress reporting (parsing files, setting up Python, applying baseline correction)
- Sample-spectrum linking - `parse_spectrum_files` command now accepts `sample_id` parameter
- `add_spectra_to_sample` method in `SampleStorage` to track spectrum ownership
- Frontend sample management UI with sidebar navigation and header editing
- Molecule dropdown selectors for Raman (DTNB, MBA, TFMBA) and Target (IgG, BSA, HER2, EpCAM, TROP2) molecules
- In-memory spectrum storage backend (`SpectrumStore`) for persisting parsed spectra
- Real-time event emission during spectrum parsing for incremental UI updates
- Event listeners in frontend store for handling spectrum parsing events
- Automatic spectrum cleanup when samples are deleted
- Enhanced file parsing validation and error handling

### Changed

- Moved file parsing code from `spectrum.rs` to `spectrum_importer.rs`
- Updated Rust structs to use `#[serde(rename_all = "camelCase")]` for JavaScript compatibility
- Baseline correction now uses streaming results via Rust channels and iterators
- Backend initiates baseline correction instead of frontend
- Updated spectrum struct and related data handling
- Split Python baseline correction into two files: core algorithm (`baseline_correction.py`) and batch processing (`batch_processor.py`)
- Spectrum import now streams results incrementally instead of returning all at once
- Moved spectrum loading logic from page component to samples store
- File drag-and-drop disabled until a sample is selected

### Fixed

- Spectrum clearing behavior for empty samples
- Spectrum-sample relationship to ensure correct merging during import

## [0.2.0] - 2025-09-01

### Added

#### Baseline Correction

- Implemented ALS (Asymmetric Least Squares) baseline correction algorithm
- Python-based implementation using scipy for accurate scientific computations
- Real-time baseline correction visualization in the UI
- Adjustable parameters for fine-tuning correction

#### Automatic Python Runtime Setup

- Runtime installation of Python environment on first launch
- Automatic download of uv package manager
- Automatic Python 3.13 installation with numpy and scipy
- No manual Python installation required by users
- Progress indicators during setup process
- Minimal app bundle size - Python downloaded only when needed

#### Enhanced Visualization

- Display both original and baseline-corrected spectra
- Visual comparison of spectra before and after correction
- Better chart formatting and presentation

### Changed

#### Build Process

- Simplified CI/CD pipeline
- Improved cross-platform build support

### Fixed

- TypeScript type errors in Svelte components (`onMount` async handling)
- Windows compatibility for Python/uv extraction process
- Cross-platform path handling for Python executables

### Developer Experience

- Added format commands for both Rust and TypeScript/Svelte code
- Improved error messages for Python setup issues
- Better documentation of the new architecture in CLAUDE.md and README.md

[0.2.0]: https://github.com/mdoverhag/raman-tools/releases/tag/v0.2.0

## [0.1.0] - 2025-08-31

### Added

#### Core Features

- Initial release of Raman Tools desktop application
- Drag-and-drop interface for bulk spectrum file uploads
- Support for .txt format spectrum files (tab-separated wavenumber/intensity data)
- Real-time spectrum visualization using D3.js
- Dark mode UI optimized for scientific data analysis

#### Technical Implementation

- Cross-platform desktop app built with Tauri v2 and SvelteKit
- High-performance Rust backend for file parsing
- Support for processing 150+ spectrum files simultaneously
- Spectrum data structure: ~1801 data points (200-2000 cm⁻¹ range)

#### Build & Distribution

- Automated CI/CD pipeline with GitHub Actions
- macOS code signing with Apple Developer ID certificate
- macOS notarization for releases (no security warnings)
- Windows MSI installer (unsigned - requires SmartScreen bypass)
- Styled DMG installer for macOS with drag-to-Applications interface

#### Development Setup

- TypeScript and Rust testing infrastructure
- Prettier and Rust fmt code formatting
- Comprehensive documentation (README.md and CLAUDE.md)
- Inter font bundled locally for consistent typography

### Known Issues

- Windows installer shows SmartScreen warning (code signing certificate pending)
- Limited to viewing single spectrum at a time
- No data export functionality yet

### Notes

This is the foundation release focused on establishing the core architecture and basic spectrum viewing capabilities. Future releases will add advanced analysis features including baseline correction, peak detection, and statistical analysis tools.

[0.1.0]: https://github.com/mdoverhag/raman-tools/releases/tag/v0.1.0
