# Changelog

All notable changes to Raman Tools will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
