use crate::python_bridge;
use crate::spectrum::{self, Spectrum};
use serde::{Deserialize, Serialize};
use tauri::{AppHandle, Emitter};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum ImportEvent {
    Progress {
        stage: String,
        current: usize,
        total: usize,
        filename: String,
    },
    SpectrumReady {
        spectrum: Spectrum,
    },
    SpectrumUpdated {
        spectrum: Spectrum,
    },
    Complete {
        spectra: Vec<Spectrum>,
    },
    Error {
        filename: String,
        error: String,
    },
}

pub async fn import_spectra(
    app: AppHandle,
    filepaths: Vec<String>,
) -> Result<Vec<Spectrum>, String> {
    let total_files = filepaths.len();
    let mut spectra = Vec::new();

    // Stage 1: Parse files
    for (index, filepath) in filepaths.iter().enumerate() {
        let filename = std::path::Path::new(filepath)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Emit parsing progress
        app.emit(
            "import:progress",
            ImportEvent::Progress {
                stage: "parsing".to_string(),
                current: index + 1,
                total: total_files,
                filename: filename.clone(),
            },
        )
        .map_err(|e| format!("Failed to emit progress event: {}", e))?;

        // Parse single file
        match spectrum::parse_files(vec![filepath.clone()]) {
            Ok(mut parsed) => {
                if let Some(spectrum) = parsed.pop() {
                    // Emit spectrum ready event immediately after parsing
                    app.emit(
                        "import:spectrum_ready",
                        ImportEvent::SpectrumReady {
                            spectrum: spectrum.clone(),
                        },
                    )
                    .map_err(|e| format!("Failed to emit spectrum ready event: {}", e))?;

                    spectra.push(spectrum);
                }
            }
            Err(e) => {
                // Emit error event but continue with other files
                app.emit(
                    "import:error",
                    ImportEvent::Error {
                        filename: filename.clone(),
                        error: e.clone(),
                    },
                )
                .map_err(|e| format!("Failed to emit error event: {}", e))?;

                eprintln!("Failed to parse {}: {}", filename, e);
            }
        }
    }

    // Stage 2: Apply baseline correction
    let total_spectra = spectra.len();

    // Baseline correction parameters (default for now)
    let denoise = true;
    let window_size = 5;
    let lambda_param = 1e7;
    let p = 0.01;
    let d = 2;

    for (index, spectrum) in spectra.iter_mut().enumerate() {
        // Emit baseline correction progress
        app.emit(
            "import:progress",
            ImportEvent::Progress {
                stage: "baseline".to_string(),
                current: index + 1,
                total: total_spectra,
                filename: spectrum.filename.clone(),
            },
        )
        .map_err(|e| format!("Failed to emit progress event: {}", e))?;

        // Convert u16 intensities to f64 for baseline correction
        let intensities_f64: Vec<f64> = spectrum.intensities.iter().map(|&i| i as f64).collect();

        match python_bridge::apply_baseline_correction(
            app.clone(),
            intensities_f64,
            denoise,
            window_size,
            lambda_param,
            p,
            d,
        )
        .await
        {
            Ok(result) => {
                spectrum.baseline = Some(result.baseline);
                spectrum.corrected = Some(result.corrected);

                // Emit spectrum updated event with baseline correction
                app.emit(
                    "import:spectrum_updated",
                    ImportEvent::SpectrumUpdated {
                        spectrum: spectrum.clone(),
                    },
                )
                .map_err(|e| format!("Failed to emit spectrum updated event: {}", e))?;
            }
            Err(e) => {
                // Emit error event but continue
                app.emit(
                    "import:error",
                    ImportEvent::Error {
                        filename: spectrum.filename.clone(),
                        error: format!("Baseline correction failed: {}", e),
                    },
                )
                .map_err(|e| format!("Failed to emit error event: {}", e))?;

                eprintln!(
                    "Failed to apply baseline correction to {}: {}",
                    spectrum.filename, e
                );
            }
        }
    }

    // Emit complete event
    app.emit(
        "import:complete",
        ImportEvent::Complete {
            spectra: spectra.clone(),
        },
    )
    .map_err(|e| format!("Failed to emit complete event: {}", e))?;

    Ok(spectra)
}
