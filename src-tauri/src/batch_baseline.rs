//! Batch baseline correction module for processing multiple spectra efficiently.
//!
//! This module communicates with a Python worker process that stays alive
//! for the entire batch, avoiding the overhead of starting Python for each spectrum.

use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc::{self, Receiver};
use std::thread;
use tauri::AppHandle;

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

/// Progress update from batch processing
#[derive(Debug, Clone)]
pub struct BatchProgress {
    pub current: usize,
    pub total: usize,
}

/// Update from batch processor - either progress or a result
#[derive(Debug)]
pub enum BatchUpdate {
    Progress(BatchProgress),
    Result(Result<SpectrumResult, String>),
}

/// Holds the Python child process and receiver for cleanup
pub struct BatchProcessor {
    child: Option<Child>,
    receiver: Receiver<BatchUpdate>,
}

impl Drop for BatchProcessor {
    fn drop(&mut self) {
        // Ensure Python process is terminated when iterator is dropped
        if let Some(mut child) = self.child.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

impl Iterator for BatchProcessor {
    type Item = BatchUpdate;

    fn next(&mut self) -> Option<Self::Item> {
        // Receive next update from channel
        match self.receiver.recv() {
            Ok(update) => Some(update),
            Err(_) => None, // Channel closed, iteration complete
        }
    }
}

/// Process multiple spectra in a single Python session, returning an iterator
///
/// This function:
/// 1. Starts a single Python process with the batch processor script
/// 2. Sends all spectra data at once
/// 3. Returns an iterator that yields results as they stream back
/// 4. Emits progress events as processing occurs
pub fn process_batch_streaming(
    app: AppHandle,
    spectra: Vec<Vec<f64>>,
    params: BaselineParams,
) -> Result<BatchProcessor, String> {
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
        stdin
            .flush()
            .map_err(|e| format!("Failed to flush stdin: {}", e))?;
    }
    child.stdin.take(); // Drop stdin to close it

    // Take stdout for processing
    let stdout = child.stdout.take().ok_or("Failed to get stdout handle")?;

    // Create channel for results
    let (sender, receiver) = mpsc::channel();

    // Spawn thread to read Python output and send to channel
    thread::spawn(move || {
        let reader = BufReader::new(stdout);

        for line in reader.lines() {
            match line {
                Ok(line) => {
                    // Parse the JSON message
                    match serde_json::from_str::<BatchMessage>(&line) {
                        Ok(message) => {
                            match message {
                                BatchMessage::Progress { current, total } => {
                                    // Send progress update through channel
                                    let _ = sender.send(BatchUpdate::Progress(BatchProgress {
                                        current,
                                        total,
                                    }));
                                }

                                BatchMessage::Result {
                                    index,
                                    baseline,
                                    corrected,
                                } => {
                                    // Send result through channel
                                    let _ = sender.send(BatchUpdate::Result(Ok(SpectrumResult {
                                        index,
                                        baseline,
                                        corrected,
                                    })));
                                }

                                BatchMessage::Error {
                                    index,
                                    error,
                                    error_type,
                                } => {
                                    // Send error through channel
                                    let _ = sender.send(BatchUpdate::Result(Err(format!(
                                        "Spectrum {} failed ({}): {}",
                                        index, error_type, error
                                    ))));
                                }

                                BatchMessage::Complete {} => {
                                    // Processing complete, close channel
                                    break;
                                }

                                BatchMessage::FatalError { error, error_type } => {
                                    // Send fatal error and close
                                    let _ = sender.send(BatchUpdate::Result(Err(format!(
                                        "Fatal Python error ({}): {}",
                                        error_type, error
                                    ))));
                                    break;
                                }
                            }
                        }
                        Err(_) => {
                            // Ignore parse errors, could be stderr output
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Failed to read output line: {}", e);
                    break;
                }
            }
        }

        // Channel will close when sender is dropped
    });

    Ok(BatchProcessor {
        child: Some(child),
        receiver,
    })
}
