# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Tauri desktop application with a SvelteKit frontend using TypeScript. The project uses Bun as the package manager and Vite as the build tool.

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
