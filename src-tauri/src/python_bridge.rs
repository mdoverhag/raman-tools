//! Python Bridge Module for Baseline Correction
//!
//! This module provides a bridge between Rust/Tauri and Python for performing
//! baseline correction on Raman spectra. It uses a runtime-installed Python with
//! numpy and scipy to execute the ALS (Asymmetric Least Squares) algorithm.
//!
//! ## Why this exists:
//! - Baseline correction algorithms like ALS require complex matrix operations
//! - Python's scipy provides optimized sparse matrix solvers that would be complex to reimplement in Rust
//! - Scientists are familiar with Python implementations and can verify/modify the algorithm
//! - Python is installed at runtime using uv for minimal bundle size
//!
//! ## Architecture:
//! - Python runtime is downloaded and installed on first run
//! - Communication via JSON through subprocess stdin/stdout
//! - Platform-specific Python paths for Windows, macOS, and Linux

use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use tauri::AppHandle;

// Import python_setup for runtime Python
use crate::python_setup;

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

/// Request structure for averaging multiple spectra
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AverageRequest {
    pub raw_spectra: Vec<Vec<f64>>,
    pub corrected_spectra: Vec<Vec<f64>>,
}

/// Response from Python containing the averaged spectra and statistics
#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AverageResponse {
    pub average_intensities: Vec<f64>,
    pub std_dev_intensities: Vec<f64>,
    pub average_corrected: Option<Vec<f64>>,
    pub std_dev_corrected: Option<Vec<f64>>,
    pub count: usize,
}

/// Enum to handle average calculation results
#[derive(Debug, Deserialize)]
#[serde(untagged)]
pub enum AverageResult {
    Success(AverageResponse),
    Error { error: String },
}

/// Get the path to Python executable from runtime-installed Python
fn get_python_path(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_python = python_setup::get_python_path(app)?;

    if !runtime_python.exists() {
        return Err(format!(
            "Python runtime not found at: {:?}. Please ensure Python is installed via the app.",
            runtime_python
        ));
    }

    Ok(runtime_python)
}

/// Get the path to the baseline correction script from runtime installation
fn get_script_path(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_script = python_setup::get_baseline_script_path(app)?;

    if !runtime_script.exists() {
        return Err(format!(
            "Baseline correction script not found at: {:?}. Please ensure Python environment is properly set up.",
            runtime_script
        ));
    }

    Ok(runtime_script)
}

/// Apply baseline correction to a spectrum using the runtime-installed Python
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

/// Calculate average and standard deviation of raw and corrected spectra
///
/// This function uses Python to efficiently calculate the mean and standard deviation
/// across multiple spectra using numpy's optimized operations.
///
/// ## Parameters:
/// - `raw_spectra`: Vector of raw intensity spectra
/// - `corrected_spectra`: Vector of baseline-corrected spectra
///
/// ## Returns:
/// - `AverageResponse` containing the averaged spectra, standard deviations, and count
/// - Error string if Python execution fails or spectra are invalid
pub async fn calculate_average(
    app: AppHandle,
    raw_spectra: Vec<Vec<f64>>,
    corrected_spectra: Vec<Vec<f64>>,
) -> Result<AverageResponse, String> {
    if raw_spectra.is_empty() {
        return Err("No spectra provided for averaging".to_string());
    }

    // Verify all raw spectra have the same length
    let expected_len = raw_spectra[0].len();
    for (i, spectrum) in raw_spectra.iter().enumerate() {
        if spectrum.len() != expected_len {
            return Err(format!(
                "Raw spectrum {} has different length ({} vs {})",
                i,
                spectrum.len(),
                expected_len
            ));
        }
    }

    // Verify corrected spectra if provided
    if !corrected_spectra.is_empty() {
        for (i, spectrum) in corrected_spectra.iter().enumerate() {
            if spectrum.len() != expected_len {
                return Err(format!(
                    "Corrected spectrum {} has different length ({} vs {})",
                    i,
                    spectrum.len(),
                    expected_len
                ));
            }
        }
    }

    // Get paths
    let python_path = get_python_path(&app)?;
    let script_path = python_setup::get_calc_averages_path(&app)?;

    // Prepare the request
    let request = AverageRequest {
        raw_spectra,
        corrected_spectra,
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
    let result: AverageResult = serde_json::from_slice(&output.stdout).map_err(|e| {
        let stdout = String::from_utf8_lossy(&output.stdout);
        format!("Failed to parse Python output: {}. Output: {}", e, stdout)
    })?;

    match result {
        AverageResult::Success(response) => Ok(response),
        AverageResult::Error { error } => Err(format!("Python error: {}", error)),
    }
}
