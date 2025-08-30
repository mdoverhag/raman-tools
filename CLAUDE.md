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

### Tauri Development
- `bun run tauri dev` - Start the Tauri app in development mode (runs both frontend and Rust backend)
- `bun run tauri build` - Build the Tauri app for production
- `cargo build` - Build the Rust backend only (run in src-tauri directory)
- `cargo test` - Run Rust tests (run in src-tauri directory)

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

### Platform-Specific
- File system access patterns differ between Windows/Mac
- Consider native file dialogs via Tauri APIs
- Ensure consistent performance across platforms
