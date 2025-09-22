# Raman Tools

A cross-platform desktop application for analyzing Raman spectroscopy data, supporting research on detecting circulating tumor cells (CTCs) in blood for breast cancer diagnosis.

## Overview

Raman Tools is designed to help scientists analyze Raman spectroscopy output files on their Windows and Mac machines. The application handles bulk processing of spectral data, visualization, and analysis of molecular interactions between Raman-active compounds and target biomolecules.

## Research Context

This tool supports PhD research focused on:

- Detecting circulating tumor cells (CTCs) in blood samples
- Identifying breast cancer biomarkers (HER2, EpCAM, TROP2)
- Analyzing interactions between Raman molecules (DTNB, MBA, TFMBA) and target proteins (IgG, BSA)
- Processing large datasets with typically 150 replicate measurements per sample

## Key Features

### Implemented Features

- **Sample Management**: Organize spectra into samples with molecular configurations
  - Create and manage samples with custom names
  - Assign Raman molecules (DTNB, MBA, TFMBA) and target molecules (IgG, BSA, HER2, EpCAM, TROP2)
  - Track spectrum count per sample
- **Bulk File Upload**: Drag-and-drop support for 150+ spectrum files with real-time progress
  - Automatic linking to selected sample
  - Progress tracking for parsing and baseline correction
- **Spectrum Parser**: Imports .txt format files with wavenumber/intensity data
- **Baseline Correction**: ALS (Asymmetric Least Squares) algorithm via integrated Python runtime
  - Single Python process for all spectra (optimized performance)
  - Streaming results with real-time UI updates
  - Automatic Python environment setup on first launch
- **Interactive Visualization**: Real-time spectrum plotting with Chart.js
  - Display raw intensities, baseline, and corrected spectra
  - Zoom and pan functionality

### In Development

- **Sample Organization**: Enhanced support for single and multiplex samples
  - Single: One Raman molecule + one target (e.g., DTNB + IgG)
  - Multiplex: Multiple molecules for detecting multiple biomarkers simultaneously
- **Experiment Tracking**: Organize samples into experiments with metadata
- **Peak Detection**: Identify and quantify characteristic peaks
- **Statistical Analysis**: Process replicate measurements
- **Batch Visualization**: Compare multiple spectra overlaid

## Technology Stack

- **Frontend**: SvelteKit + TypeScript (Svelte 5)
- **Backend**: Rust + Tauri
- **Build Tool**: Vite
- **Package Manager**: Bun
- **Platform**: Cross-platform desktop (Windows, macOS, Linux)

## Installation

### Download Pre-built Binaries

Download the latest release from the [GitHub Releases](https://github.com/yourusername/raman-tools/releases) page:

- **Windows**: `Raman Tools_[version]_x64_en-US.msi` - MSI installer for Windows x86_64
- **macOS**: `Raman Tools_[version]_aarch64.dmg` - DMG installer for Apple Silicon Macs

#### Installation Notes

**Windows**: The installer is not code signed (certificates cost $200-600/year). Windows Defender SmartScreen will show a warning. To install:

1. Click "More info" on the SmartScreen warning
2. Click "Run anyway"

The installer is built automatically by GitHub Actions from public source code.

> **Note**: Microsoft's [Azure Trusted Signing](https://azure.microsoft.com/en-us/products/trusted-signing) service may offer affordable signing for individual developers in the future. Currently restricted to organizations with 3+ year history.

**macOS**: The app is signed with an Apple Developer ID certificate and notarized by Apple. It will open without any security warnings.

**First Launch**: The application will automatically download and configure Python for baseline correction on first run. This requires an internet connection and takes about 1-2 minutes.

## Development

### Prerequisites

- [Rust](https://www.rust-lang.org/tools/install)
- [Bun](https://bun.sh)
- [Node.js](https://nodejs.org) (for tooling compatibility)

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/raman-tools.git
cd raman-tools

# Install dependencies
bun install

# Start development server
bun run tauri dev
```

**Note**: On first launch, the application will automatically download and set up Python with numpy and scipy for baseline correction. This is a one-time setup that requires an internet connection.

### Building

```bash
# Build for production
bun run tauri build

# Build for specific target
bun run tauri build --target x86_64-pc-windows-msvc  # Windows
bun run tauri build --target aarch64-apple-darwin    # macOS Apple Silicon
```

### CI/CD

The project uses GitHub Actions for continuous integration and deployment:

- **Tests**: Run on every PR and push to master (TypeScript checks, Rust tests, formatting)
- **Builds**: Create installers for Windows and macOS on every push using `tauri-action`
- **Releases**: Automatically create GitHub releases with installers when pushing version tags

To create a new release:

```bash
git tag v0.5.0
git push origin v0.5.0
```

## Data Format

The application expects spectrum files in .txt format with:

- Wavenumber range: typically 200-2000 cm⁻¹
- Step size: typically 1.0 cm⁻¹
- Two columns: wavenumber and intensity
- Files named sequentially (e.g., Spectrum_001.txt, Spectrum_002.txt)

## Roadmap

### Phase 1: Core Infrastructure ✅

- [x] Basic Tauri application setup
- [x] File system integration for bulk uploads
- [x] Local data storage (in-memory with planned SQLite migration)

### Phase 2: Data Import & Visualization ✅

- [x] Spectrum file parser
- [x] Interactive plotting with Chart.js
- [x] Basic sample management UI
- [x] Real-time progress tracking for imports

### Phase 3: Analysis Features (In Progress)

- [x] Baseline correction with ALS algorithm
- [x] Python integration with automatic environment setup
- [ ] Peak detection and quantification
- [ ] Export processed data for Python/R analysis

### Phase 4: Advanced Features

- [ ] Statistical analysis of replicates
- [ ] Machine learning integration for pattern recognition
- [ ] Report generation

## Contributing

This project is currently in early development. Contributions and suggestions are welcome!

## License

MIT
