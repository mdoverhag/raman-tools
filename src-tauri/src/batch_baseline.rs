//! Batch baseline correction module for processing multiple spectra efficiently.
//!
//! This module communicates with a Python worker process that stays alive
//! for the entire batch, avoiding the overhead of starting Python for each spectrum.

use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::process::{Command, Stdio};
use tauri::{AppHandle, Emitter};

use crate::python_setup;

/// Parameters for baseline correction
#[derive(Debug, Clone, Serialize)]
pub struct BaselineParams {
    pub denoise: bool,
    pub window_size: i32,
    pub lambda_param: f64,
    pub p: f64,
    pub d: i32,
}

impl Default for BaselineParams {
    fn default() -> Self {
        Self {
            denoise: true,
            window_size: 5,
            lambda_param: 1e7,
            p: 0.01,
            d: 2,
        }
    }
}

/// Input spectrum for batch processing
#[derive(Debug, Serialize)]
struct InputSpectrum {
    index: usize,
    intensities: Vec<f64>,
}

/// Batch request sent to Python
#[derive(Debug, Serialize)]
struct BatchRequest {
    spectra: Vec<InputSpectrum>,
    params: BaselineParams,
}

/// Message types from Python
#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
enum BatchMessage {
    #[serde(rename = "progress")]
    Progress { current: usize, total: usize },

    #[serde(rename = "result")]
    Result {
        index: usize,
        baseline: Vec<f64>,
        corrected: Vec<f64>,
    },

    #[serde(rename = "error")]
    Error {
        index: usize,
        error: String,
        error_type: String,
    },

    #[serde(rename = "complete")]
    Complete {},

    #[serde(rename = "fatal_error")]
    FatalError { error: String, error_type: String },
}

/// Result for a single spectrum
#[derive(Debug, Clone)]
pub struct SpectrumResult {
    pub index: usize,
    pub baseline: Vec<f64>,
    pub corrected: Vec<f64>,
}

/// Process multiple spectra in a single Python session
///
/// This function:
/// 1. Starts a single Python process with the batch processor script
/// 2. Sends all spectra data at once
/// 3. Reads results as they stream back, emitting progress events
/// 4. Returns all results when complete
pub async fn process_batch(
    app: AppHandle,
    spectra: Vec<Vec<f64>>,
    params: BaselineParams,
    progress_event: &str,
) -> Result<Vec<SpectrumResult>, String> {
    // Get Python and script paths
    let python_path = python_setup::get_python_path(&app)?;
    let batch_script_path = python_setup::get_batch_processor_path(&app)?;

    if !batch_script_path.exists() {
        return Err(format!(
            "Batch processor script not found at: {:?}",
            batch_script_path
        ));
    }

    // Prepare batch request
    let input_spectra: Vec<InputSpectrum> = spectra
        .into_iter()
        .enumerate()
        .map(|(index, intensities)| InputSpectrum { index, intensities })
        .collect();

    let request = BatchRequest {
        spectra: input_spectra,
        params,
    };

    let request_json = serde_json::to_string(&request)
        .map_err(|e| format!("Failed to serialize batch request: {}", e))?;

    eprintln!(
        "Sending batch request with {} spectra",
        request.spectra.len()
    );
    eprintln!(
        "First spectrum has {} intensities",
        request
            .spectra
            .first()
            .map(|s| s.intensities.len())
            .unwrap_or(0)
    );

    // Start Python process
    let mut child = Command::new(&python_path)
        .arg(&batch_script_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start Python batch processor: {}", e))?;

    // Send the batch request
    {
        let stdin = child.stdin.as_mut().ok_or("Failed to get stdin handle")?;
        stdin
            .write_all(request_json.as_bytes())
            .map_err(|e| format!("Failed to write batch request: {}", e))?;
        // Close stdin to signal we're done sending
        stdin
            .flush()
            .map_err(|e| format!("Failed to flush stdin: {}", e))?;
    }
    child.stdin.take(); // Drop stdin to close it

    // Read streaming results
    let stdout = child.stdout.take().ok_or("Failed to get stdout handle")?;
    let reader = BufReader::new(stdout);
    let mut results = Vec::new();
    let mut errors = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Failed to read output line: {}", e))?;

        // Parse the JSON message
        let message: BatchMessage = serde_json::from_str(&line)
            .map_err(|e| format!("Failed to parse message: {}. Line: {}", e, line))?;

        match message {
            BatchMessage::Progress { current, total } => {
                // Emit progress event if event name provided
                if !progress_event.is_empty() {
                    app.emit(
                        progress_event,
                        serde_json::json!({
                            "stage": "baseline",
                            "current": current,
                            "total": total,
                            "filename": format!("Processing {} of {}", current, total),
                        }),
                    )
                    .ok(); // Ignore emission errors
                }
            }

            BatchMessage::Result {
                index,
                baseline,
                corrected,
            } => {
                results.push(SpectrumResult {
                    index,
                    baseline,
                    corrected,
                });
            }

            BatchMessage::Error {
                index,
                error,
                error_type,
            } => {
                errors.push(format!(
                    "Spectrum {} failed ({}): {}",
                    index, error_type, error
                ));
            }

            BatchMessage::Complete {} => {
                break;
            }

            BatchMessage::FatalError { error, error_type } => {
                return Err(format!("Fatal Python error ({}): {}", error_type, error));
            }
        }
    }

    // Wait for process to complete
    let output = child
        .wait_with_output()
        .map_err(|e| format!("Failed to wait for Python process: {}", e))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("Python process failed: {}", stderr));
    }

    // Report any individual errors
    if !errors.is_empty() {
        eprintln!("Some spectra failed processing:\n{}", errors.join("\n"));
    }

    // Sort results by index to maintain order
    results.sort_by_key(|r| r.index);

    Ok(results)
}
