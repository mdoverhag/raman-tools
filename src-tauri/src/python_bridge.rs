//! Python Bridge Module for Baseline Correction
//!
//! This module provides a bridge between Rust/Tauri and Python for performing
//! baseline correction on Raman spectra. It uses a bundled Python runtime with
//! numpy and scipy to execute the ALS (Asymmetric Least Squares) algorithm.
//!
//! ## Why this exists:
//! - Baseline correction algorithms like ALS require complex matrix operations
//! - Python's scipy provides optimized sparse matrix solvers that would be complex to reimplement in Rust
//! - Scientists are familiar with Python implementations and can verify/modify the algorithm
//! - Bundling Python ensures users don't need to install Python or manage dependencies
//!
//! ## Architecture:
//! - Python runtime is bundled with the app (no user installation required)
//! - Communication via JSON through subprocess stdin/stdout
//! - Platform-specific Python paths for Windows, macOS, and Linux

use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use tauri::{AppHandle, Manager};

/// Request structure sent to Python baseline correction script
/// Contains the spectrum data and all parameters for the ALS algorithm
#[derive(Debug, Serialize)]
pub struct BaselineRequest {
    pub spectrum: Vec<f64>,
    pub denoise: bool,
    pub window_size: i32,
    pub lambda_param: f64,
    pub p: f64,
    pub d: i32,
}

/// Response from Python containing the corrected spectrum and extracted baseline
/// The corrected spectrum has the baseline removed (corrected = original - baseline)
#[derive(Debug, Deserialize, Serialize)]
pub struct BaselineResponse {
    pub corrected: Vec<f64>,
    pub baseline: Vec<f64>,
    pub denoised: Option<Vec<f64>>,
}

/// Enum to handle both successful results and Python errors
#[derive(Debug, Deserialize)]
#[serde(untagged)]
pub enum PythonResult {
    Success(BaselineResponse),
    Error { error: String, r#type: String },
}

/// Get the path to the bundled Python executable
///
/// This function locates the Python runtime that was bundled with the application
/// during the build process. The location varies by platform:
/// - macOS: Resources/python/macos/python
/// - Windows: Resources/python/windows/Scripts/python.exe
/// - Linux: Resources/python/linux/bin/python3
///
/// In development mode, it falls back to the development bundle location
fn get_python_path(app: &AppHandle) -> Result<PathBuf, String> {
    let resource_dir = app
        .path()
        .resource_dir()
        .map_err(|e| format!("Failed to get resource directory: {}", e))?;

    // Platform-specific Python location
    let python_path = if cfg!(target_os = "macos") {
        resource_dir.join("python").join("macos").join("python")
    } else if cfg!(target_os = "windows") {
        resource_dir
            .join("python")
            .join("windows")
            .join("Scripts")
            .join("python.exe")
    } else {
        return Err("Unsupported platform".to_string());
    };

    if !python_path.exists() {
        // Development fallback - currently only supports macOS
        if cfg!(debug_assertions) {
            // In development, we only support macOS
            // This fallback is hardcoded to macOS paths and won't work on Windows
            let dev_python = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("resources")
                .join("python")
                .join("macos")
                .join("python");
            if dev_python.exists() {
                return Ok(dev_python);
            }
            // No fallback to system Python - bundle must exist
            return Err("Python bundle not found. Run build-python.sh to create it.".to_string());
        }
        return Err(format!("Python runtime not found at: {:?}", python_path));
    }

    Ok(python_path)
}

/// Get the path to the baseline correction script
///
/// Locates the baseline_correction.py script that implements the ALS algorithm.
/// The script is bundled alongside the Python runtime in platform-specific locations.
/// Falls back to development paths when running in debug mode.
fn get_script_path(app: &AppHandle) -> Result<PathBuf, String> {
    let resource_dir = app
        .path()
        .resource_dir()
        .map_err(|e| format!("Failed to get resource directory: {}", e))?;

    // Platform-specific script location
    let script_path = if cfg!(target_os = "macos") {
        resource_dir
            .join("python")
            .join("macos")
            .join("baseline_correction.py")
    } else if cfg!(target_os = "windows") {
        resource_dir
            .join("python")
            .join("windows")
            .join("baseline_correction.py")
    } else {
        return Err("Unsupported platform".to_string());
    };

    // Development fallback - currently only supports macOS
    if !script_path.exists() && cfg!(debug_assertions) {
        // In development, we only support macOS for now
        // Try the bundled location first
        let bundled_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("resources")
            .join("python")
            .join("macos")
            .join("baseline_correction.py");
        if bundled_path.exists() {
            return Ok(bundled_path);
        }

        // Then try the source location
        let dev_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("python")
            .join("baseline_correction.py");
        if dev_path.exists() {
            return Ok(dev_path);
        }
    }

    if !script_path.exists() {
        return Err(format!(
            "Baseline correction script not found at: {:?}",
            script_path
        ));
    }

    Ok(script_path)
}

/// Apply baseline correction to a spectrum using the bundled Python runtime
///
/// This is the main entry point for baseline correction. It:
/// 1. Locates the bundled Python runtime and baseline correction script
/// 2. Serializes the spectrum and parameters to JSON
/// 3. Spawns a Python subprocess and sends data via stdin
/// 4. Receives the corrected spectrum and baseline via stdout
/// 5. Handles errors from both the subprocess and Python script
///
/// ## Parameters:
/// - `spectrum`: Raw intensity values from the Raman spectrum
/// - `denoise`: Whether to apply moving average denoising before baseline correction
/// - `window_size`: Size of the moving average window (if denoising)
/// - `lambda_param`: ALS smoothness parameter (larger = smoother baseline, typical: 1e7)
/// - `p`: ALS asymmetry parameter (typical: 0.01 for Raman spectra)
/// - `d`: Order of differences for the penalty matrix (typically 2)
///
/// ## Returns:
/// - `BaselineResponse` containing corrected spectrum and baseline
/// - Error string if Python execution fails
pub async fn apply_baseline_correction(
    app: AppHandle,
    spectrum: Vec<f64>,
    denoise: bool,
    window_size: i32,
    lambda_param: f64,
    p: f64,
    d: i32,
) -> Result<BaselineResponse, String> {
    // Get paths
    let python_path = get_python_path(&app)?;
    let script_path = get_script_path(&app)?;

    // Prepare the request
    let request = BaselineRequest {
        spectrum,
        denoise,
        window_size,
        lambda_param,
        p,
        d,
    };

    let request_json = serde_json::to_string(&request)
        .map_err(|e| format!("Failed to serialize request: {}", e))?;

    // Run Python subprocess
    let mut child = Command::new(&python_path)
        .arg(&script_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| {
            format!(
                "Failed to start Python process: {}. Python path: {:?}",
                e, python_path
            )
        })?;

    // Send input
    {
        let stdin = child.stdin.as_mut().ok_or("Failed to get stdin handle")?;
        stdin
            .write_all(request_json.as_bytes())
            .map_err(|e| format!("Failed to write to Python process: {}", e))?;
    }

    // Wait for output
    let output = child
        .wait_with_output()
        .map_err(|e| format!("Failed to wait for Python process: {}", e))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("Python process failed: {}", stderr));
    }

    // Parse the response
    let result: PythonResult = serde_json::from_slice(&output.stdout).map_err(|e| {
        let stdout = String::from_utf8_lossy(&output.stdout);
        format!("Failed to parse Python output: {}. Output: {}", e, stdout)
    })?;

    match result {
        PythonResult::Success(response) => Ok(response),
        PythonResult::Error { error, r#type } => {
            Err(format!("Python error ({}): {}", r#type, error))
        }
    }
}

/// Check if Python runtime is available
///
/// Verifies that both the Python executable and baseline correction script
/// are present and accessible. Used by the frontend to check if baseline
/// correction functionality is available before attempting to use it.
pub fn check_python_availability(app: &AppHandle) -> Result<bool, String> {
    let python_path = get_python_path(app)?;
    let script_path = get_script_path(app)?;

    Ok(python_path.exists() && script_path.exists())
}

/// Get Python runtime information for diagnostics
///
/// Returns version and path information about the bundled Python runtime.
/// Useful for debugging issues and verifying the correct Python is being used.
pub fn get_python_info(app: &AppHandle) -> Result<String, String> {
    let python_path = get_python_path(app)?;

    let output = Command::new(&python_path)
        .arg("--version")
        .output()
        .map_err(|e| format!("Failed to get Python version: {}", e))?;

    let version = String::from_utf8_lossy(&output.stdout);
    Ok(format!(
        "Python path: {:?}\nVersion: {}",
        python_path, version
    ))
}
