# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Changelog Style

When updating CHANGELOG.md:

- Use technical, straightforward language - describe what actually changed in the code
- Avoid marketing-style language or embellishments
- List each change only once in the most appropriate section
- Focus on concrete changes: new files, moved code, refactored modules, etc.
- Don't duplicate items across multiple sections

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
- **Important**: Rust structs use snake_case but must serialize to camelCase for JavaScript using `#[serde(rename_all = "camelCase")]`

### Python Integration for Baseline Correction

- **Purpose**: Leverages Python's scipy for complex matrix operations (ALS algorithm)
- **Architecture**: Runtime-installed Python with numpy/scipy via uv package manager
- **Communication**: JSON via subprocess stdin/stdout between Rust and Python
- **Installation**: Automatic on first launch - downloads uv, installs Python 3.13, and sets up virtual environment
- **Platform Support**: Full support for Windows, macOS, and Linux
- **Location**: Python runtime stored in app data directory (`~/Library/Application Support/com.mikaeldoverhag.raman-tools/runtime/` on macOS)
- **Source Files**:
  - `baseline_correction.py`: Pure algorithm implementation (ALS and denoising)
  - `batch_processor.py`: Handles batch processing with streaming JSON output
  - `requirements.txt`: Python dependencies (numpy, scipy)
- **File Syncing**: Python files are automatically synced to runtime directory on app startup via `sync_python_files()`

### Batch Processing Architecture

- **Performance Optimization**: Single Python process handles all spectra (avoids ~150 process starts)
- **Streaming Results**: Uses Rust channels and iterators to stream results as they complete
- **Real-time Updates**: UI updates immediately as each spectrum is processed
- **Incremental Processing**: Spectra are added to backend storage as they're parsed, not in batch
- **Progress Tracking**: Three-stage progress reporting:
  1. Parsing files (shows filename, validates data)
  2. Preparing baseline correction (Python bytecode compilation on first run)
  3. Applying baseline correction (shows filename)
- **Module Separation**:
  - `batch_baseline.rs`: Manages Python process, provides iterator interface with `BatchUpdate` enum
  - `spectrum_importer.rs`: Orchestrates import flow, validates data, and emits UI events
  - `spectrum.rs`: Manages spectrum storage with add/get/delete operations
  - Clean separation between data processing and UI event emission
- **Data Validation**: Enhanced error handling and validation during file parsing

## Key Configuration Files

- **Frontend Build**: Outputs to `/build/` directory (configured in tauri.conf.json)
- **Tauri Config**: Uses Bun commands for development and build processes
- **Package Manager**: Uses Bun (not npm/yarn/pnpm)

## Sample Management System

- **Storage**: In-memory HashMap in Rust backend (`src-tauri/src/samples.rs`)
- **Spectrum Storage**: In-memory HashMap in Rust backend (`src-tauri/src/spectrum.rs`)
  - `SpectrumStore` manages spectrum persistence
  - Automatic cleanup when samples are deleted
- **Sample Model**: Each sample has:
  - ID (UUID)
  - Name
  - Raman molecules (DTNB, MBA, TFMBA)
  - Target molecules (IgG, BSA, HER2, EpCAM, TROP2)
  - Spectrum IDs (links to uploaded spectra)
- **UI Patterns**:
  - Left sidebar shows sample list (read-only)
  - Double-click sample name in header to edit
  - Click molecule pills to open dropdown selector
  - Spectra automatically link to selected sample on upload
  - File drag-and-drop disabled until sample is selected
- **Frontend State**: Svelte store (`samples.svelte.ts`) caches data and manages UI state
  - Handles spectrum loading and merging
  - Listens to real-time parsing events via Tauri event system
- **Key Flow**: Select sample → Drop files → Files linked to sample → Spectrum count updates incrementally
- **Event System**:
  - `spectrum-parsed` events emitted during import for real-time UI updates
  - Frontend store subscribes to events and updates state incrementally

## Domain Context

### Scientific Background

- **Purpose**: Analyze Raman spectroscopy data for detecting cancer biomarkers
- **Key Molecules**:
  - Raman-active compounds: DTNB, MBA, TFMBA
  - Target proteins/biomarkers: IgG, BSA, HER2, EpCAM, TROP2
- **Data Structure**: Samples typically contain 150 replicate spectrum measurements
- **File Format**: .txt files with wavenumber (200-2000 cm⁻¹) and intensity columns

### Data Model (Implemented)

- **Samples**: Container for organizing spectra
  - Has name and molecule configuration
  - Tracks which spectra belong to it
  - Can have multiple Raman and Target molecules
- **Spectra**: Individual measurements (~1801 data points each)
  - Linked to a sample via sample_id
  - Contains wavenumber/intensity arrays
  - Has baseline correction applied automatically
- **Future**: Experiments (collections of samples) not yet implemented

## Implementation Status

### Completed Features

- ✅ Basic Tauri application with dark theme UI
- ✅ Drag-and-drop bulk file upload (150+ files)
- ✅ Spectrum file parser for .txt format
- ✅ Real-time plotting with Chart.js
- ✅ Sample management system (CRUD operations)
- ✅ Baseline correction via Python integration (ALS algorithm)
- ✅ Batch processing with progress updates
- ✅ Linking spectra to samples

### In Progress

- 🔄 Data persistence (currently in-memory only)
- 🔄 Export functionality

### Future Features

- Statistical analysis of replicate measurements
- Peak detection and quantification
- Experiment management (grouping samples)
- Report generation
- SQLite or file-based persistence

## CI/CD and Releases

### GitHub Actions Workflow

The project has automated CI/CD configured in `.github/workflows/build.yml`:

1. **Test Job**: Runs on all PRs and pushes to master
   - Installs Linux dependencies (libgtk-3-dev, libwebkit2gtk-4.1-dev, etc.)
   - Runs TypeScript type checking (`bun run check`)
   - Runs Rust tests (`cargo test`)
   - Checks Rust formatting (`cargo fmt --check`)
   - Runs Rust linter (`cargo clippy`)

2. **Build Job**: Runs after tests pass
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
