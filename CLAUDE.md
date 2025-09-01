# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Raman Tools is a cross-platform desktop application for analyzing Raman spectroscopy data, built with Tauri and SvelteKit. The application supports scientific research on detecting circulating tumor cells (CTCs) in blood for breast cancer diagnosis.

## Development Commands

### Frontend Development

- `bun run dev` - Start the SvelteKit development server (runs on http://localhost:1420)
- `bun run build` - Build the frontend for production
- `bun run preview` - Preview the production build
- `bun run check` - Run svelte-kit sync and type checking
- `bun run check:watch` - Run type checking in watch mode
- `bun run format` - Format all JS/TS/Svelte files with Prettier

### Tauri Development

- `bun run tauri dev` - Start the Tauri app in development mode (runs both frontend and Rust backend)
- `bun run tauri build` - Build the Tauri app for production
- `cargo build` - Build the Rust backend only (run in src-tauri directory)
- `cargo test` - Run Rust tests (run in src-tauri directory)
- `cargo fmt` - Format Rust code (run in src-tauri directory)

## Code Formatting

**IMPORTANT**: After making code changes, always run:

1. `bun run format` - Formats JavaScript/TypeScript/Svelte files (also formats YAML files)
2. `cd src-tauri && cargo fmt` - Formats Rust files

These commands ensure consistent code style across the project.

## Architecture

### Frontend (SvelteKit + TypeScript)

- **Location**: `/src/` directory
- **Framework**: SvelteKit with Svelte 5 using runes ($state, $derived, etc.)
- **Build Tool**: Vite
- **Routing**: File-based routing in `/src/routes/`
- **Static Adapter**: Uses adapter-static for SPA mode (required for Tauri)

### Backend (Rust + Tauri)

- **Location**: `/src-tauri/` directory
- **Entry Point**: `src-tauri/src/main.rs` calls `raman_tools_lib::run()`
- **Commands**: Defined in `src-tauri/src/lib.rs` with `#[tauri::command]` attribute
- **Configuration**: `src-tauri/tauri.conf.json` defines app settings, window configuration, and build settings

### Frontend-Backend Communication

- Frontend uses `@tauri-apps/api/core` to invoke Rust commands
- Example: `await invoke("greet", { name })` calls the `greet` function in Rust
- Commands must be registered in `tauri::generate_handler![]` in `lib.rs`

### Python Integration for Baseline Correction

- **Purpose**: Leverages Python's scipy for complex matrix operations (ALS algorithm)
- **Architecture**: Bundled Python runtime with numpy/scipy, no user installation required
- **Communication**: JSON via subprocess stdin/stdout between Rust and Python
- **Build Scripts**:
  - `src-tauri/build-python.sh` - Creates Python bundle for macOS/Linux
  - `src-tauri/build-python.ps1` - Creates Python bundle for Windows
- **Platform Support**: Full support for Windows, macOS, and Linux (Linux mainly for CI)
- **Location**: Python bundles stored in `src-tauri/resources/python/{platform}/`
- **Development**: Requires running build script before `bun run tauri dev`

## Key Configuration Files

- **Frontend Build**: Outputs to `/build/` directory (configured in tauri.conf.json)
- **Tauri Config**: Uses Bun commands for development and build processes
- **Package Manager**: Uses Bun (not npm/yarn/pnpm)

## Domain Context

### Scientific Background

- **Purpose**: Analyze Raman spectroscopy data for detecting cancer biomarkers
- **Key Molecules**:
  - Raman-active compounds: DTNB, MBA, TFMBA
  - Target proteins/biomarkers: IgG, BSA, HER2, EpCAM, TROP2
- **Data Structure**: Samples typically contain 150 replicate spectrum measurements
- **File Format**: .txt files with wavenumber (200-2000 cm⁻¹) and intensity columns

### Data Model Concepts

- **Experiments**: Research projects containing multiple samples
- **Samples**: Collections of spectra, can be:
  - Single: One Raman molecule + one target
  - Multiplex: Multiple molecules for simultaneous detection
- **Spectra**: Individual measurements with wavenumber/intensity arrays (~1801 data points each)

## Implementation Roadmap

### Phase 1: Core Infrastructure (Current)

- Set up basic Tauri application structure
- Implement file system access for bulk file uploads
- Design local data storage (consider SQLite for portability)

### Phase 2: Data Import & Visualization

- Build spectrum file parser for .txt format
- Create UI for drag-and-drop bulk file upload (150+ files)
- Implement basic plotting using a JavaScript charting library
- Design sample and experiment management interface

### Phase 3: Analysis Features

- Implement baseline correction algorithms (ALS)
- Add peak detection and quantification
- Create data export functionality for Python/R workflows

### Phase 4: Advanced Features

- Statistical analysis of replicate measurements
- Integration with scientific Python libraries (via Tauri commands)
- Report generation and data export

## CI/CD and Releases

### GitHub Actions Workflow

The project has automated CI/CD configured in `.github/workflows/build.yml`:

1. **Test Job**: Runs on all PRs and pushes to master
   - Installs Linux dependencies (libgtk-3-dev, libwebkit2gtk-4.1-dev, etc.)
   - Installs uv and builds Python bundle for Linux
   - Runs TypeScript type checking (`bun run check`)
   - Runs Rust tests (`cargo test`)
   - Checks Rust formatting (`cargo fmt --check`)
   - Runs Rust linter (`cargo clippy`)

2. **Build Job**: Runs after tests pass
   - Installs uv and builds Python bundle for target platform
   - Uses `tauri-apps/tauri-action` for consistent builds
   - Builds Windows MSI installer (x86_64)
   - Builds macOS DMG installer (Apple Silicon/aarch64) with styled installer window
   - Handles code signing automatically via environment variables
   - Uploads artifacts to GitHub Actions

3. **Release Job**: Only runs on version tags (v\*)
   - Downloads build artifacts
   - Creates GitHub release with installers attached

### Creating a Release

To create a new release:

1. Update version in `src-tauri/tauri.conf.json`
2. Commit changes
3. Create and push a version tag:
   ```bash
   git tag v0.1.0
   git push origin v0.1.0
   ```
4. GitHub Actions will automatically build and create a release

### Build Targets

- **Windows**: x86_64-pc-windows-msvc (creates .msi installer)
- **macOS**: aarch64-apple-darwin (creates .dmg for Apple Silicon)

### Code Signing

#### macOS (Implemented)

- **Signing**: All builds are signed with Apple Developer ID certificate
- **Notarization**: Release builds (version tags) are notarized by Apple
- **Required Secrets**:
  - `MACOS_CERTIFICATE`: Base64-encoded .p12 certificate
  - `MACOS_CERTIFICATE_PWD`: Certificate password
  - `APPLE_SIGNING_ID`: Full signing identity string
  - `APPLE_TEAM_ID`: Team ID (e.g., 3C3KW7PU5V)
  - `APPLE_ID`: Apple ID email (for notarization)
  - `APPLE_APP_PASSWORD`: App-specific password (for notarization)

#### Windows (Not Implemented)

- Currently unsigned due to certificate costs ($200-600/year)
- Users will see SmartScreen warning and need to click "Run anyway"
- Future options:
  - **Azure Trusted Signing**: Currently restricted to US/Canada orgs with 3+ year history. Individual developer access is paused. Check [Azure Trusted Signing updates](https://azure.microsoft.com/en-us/products/trusted-signing) and the [Tech Community blog](https://techcommunity.microsoft.com/t5/microsoft-security-blog/bg-p/Microsoft-Security-Blog) for when individual access reopens
  - **SignPath.io**: Free for open source projects
  - **Traditional certificates**: $200-600/year from DigiCert, Sectigo, etc.

## Technical Considerations

### File Processing

- Must handle bulk uploads of 150+ files efficiently
- Parse files with ~1801 data points each
- Consider using Rust for performance-critical parsing

### Data Storage

- Local-first approach for offline use
- Consider SQLite for cross-platform compatibility
- Store spectra as arrays for efficient retrieval

### UI/UX Design

- Scientists expect desktop application behavior (file dialogs, drag-and-drop)
- Need responsive charts for large datasets
- Batch operations for processing multiple samples
- **Tailwind UI Plus Available**: User has access to Tailwind UI Plus components. Ask for specific components when needed for better UI patterns (e.g., "Could you provide the Tailwind UI component for [specific pattern]?")

### Platform-Specific

- File system access patterns differ between Windows/Mac
- Consider native file dialogs via Tauri APIs
- Ensure consistent performance across platforms
