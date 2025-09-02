use crate::batch_baseline::{process_batch_streaming, BaselineParams, BatchUpdate};
use crate::spectrum::Spectrum;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;
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

// Parse a single spectrum file
fn parse_file(filepath: &str) -> Result<Spectrum, String> {
    let path = Path::new(filepath);
    let filename = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown")
        .to_string();

    let content = fs::read_to_string(filepath)
        .map_err(|e| format!("Failed to read file {}: {}", filepath, e))?;

    let mut wavenumbers: Vec<f32> = Vec::new();
    let mut intensities: Vec<u16> = Vec::new();

    for line in content.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() == 2 {
            if let (Ok(wn), Ok(intensity)) = (parts[0].parse::<f32>(), parts[1].parse::<f32>()) {
                wavenumbers.push(wn);
                intensities.push(intensity as u16);
            }
        }
    }

    if wavenumbers.is_empty() {
        return Err(format!("No valid spectrum data found in {}", filename));
    }

    // Extract wavenumber parameters
    let wavenumber_start = wavenumbers[0] as u16;
    let wavenumber_end = wavenumbers.last().unwrap().round() as u16;
    let wavenumber_step = if wavenumbers.len() > 1 {
        (wavenumbers[1] - wavenumbers[0]).round() as u16
    } else {
        1
    };

    Ok(Spectrum::new(
        filename,
        wavenumber_start,
        wavenumber_end,
        wavenumber_step,
        intensities,
    ))
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
        match parse_file(filepath) {
            Ok(spectrum) => {
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

    // Stage 2: Apply baseline correction using batch processing
    if !spectra.is_empty() {
        // Emit event that we're preparing Python
        app.emit(
            "import:progress",
            ImportEvent::Progress {
                stage: "preparing".to_string(),
                current: 0,
                total: 1,
                filename: "".to_string(),
            },
        )
        .map_err(|e| format!("Failed to emit progress event: {}", e))?;

        // Prepare intensities for batch processing
        let intensities_batch: Vec<Vec<f64>> = spectra
            .iter()
            .map(|s| s.intensities.iter().map(|&i| i as f64).collect())
            .collect();

        // Use default baseline parameters
        let params = BaselineParams::default();

        // Process spectra with streaming iterator
        match process_batch_streaming(app.clone(), intensities_batch, params) {
            Ok(processor) => {
                // Process updates as they stream in
                for update in processor {
                    match update {
                        BatchUpdate::Progress(progress) => {
                            // Get the filename from the spectrum being processed
                            let filename =
                                if progress.current > 0 && progress.current <= spectra.len() {
                                    spectra[progress.current - 1].filename.clone()
                                } else {
                                    format!("Processing {} of {}", progress.current, progress.total)
                                };

                            // Emit progress event for baseline correction
                            app.emit(
                                "import:progress",
                                ImportEvent::Progress {
                                    stage: "baseline".to_string(),
                                    current: progress.current,
                                    total: progress.total,
                                    filename,
                                },
                            )
                            .map_err(|e| format!("Failed to emit progress event: {}", e))?;
                        }
                        BatchUpdate::Result(Ok(spectrum_result)) => {
                            if spectrum_result.index < spectra.len() {
                                let spectrum = &mut spectra[spectrum_result.index];
                                spectrum.baseline = Some(spectrum_result.baseline);
                                spectrum.corrected = Some(spectrum_result.corrected);

                                // Emit spectrum updated event immediately
                                app.emit(
                                    "import:spectrum_updated",
                                    ImportEvent::SpectrumUpdated {
                                        spectrum: spectrum.clone(),
                                    },
                                )
                                .map_err(|e| {
                                    format!("Failed to emit spectrum updated event: {}", e)
                                })?;
                            }
                        }
                        BatchUpdate::Result(Err(e)) => {
                            // Log individual spectrum error but continue processing
                            app.emit(
                                "import:error",
                                ImportEvent::Error {
                                    filename: "spectrum".to_string(),
                                    error: e,
                                },
                            )
                            .map_err(|e| format!("Failed to emit error event: {}", e))?;
                        }
                    }
                }
            }
            Err(e) => {
                // Emit error for batch processing setup failure
                app.emit(
                    "import:error",
                    ImportEvent::Error {
                        filename: "batch_processing".to_string(),
                        error: format!("Failed to start batch baseline correction: {}", e),
                    },
                )
                .map_err(|e| format!("Failed to emit error event: {}", e))?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn create_test_spectrum_file(dir: &TempDir, filename: &str, content: &str) -> PathBuf {
        let file_path = dir.path().join(filename);
        fs::write(&file_path, content).unwrap();
        file_path
    }

    #[test]
    fn test_parse_valid_spectrum_file() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\n201\t143\n202\t143\n203\t144";
        let file_path = create_test_spectrum_file(&temp_dir, "test_spectrum.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_ok());
        let spectrum = result.unwrap();
        assert_eq!(spectrum.filename, "test_spectrum.txt");
        assert_eq!(spectrum.wavenumber_start, 200);
        assert_eq!(spectrum.wavenumber_end, 203);
        assert_eq!(spectrum.wavenumber_step, 1);
        assert_eq!(spectrum.intensities, vec![145, 143, 143, 144]);
    }

    #[test]
    fn test_parse_empty_file() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = create_test_spectrum_file(&temp_dir, "empty.txt", "");

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        assert!(result.unwrap_err().contains("No valid spectrum data"));
    }

    #[test]
    fn test_parse_malformed_lines() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\ninvalid line\n202\t143\n203\tnot_a_number";
        let file_path = create_test_spectrum_file(&temp_dir, "malformed.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_ok());
        let spectrum = result.unwrap();
        assert_eq!(spectrum.intensities.len(), 2); // Only valid lines
    }

    #[test]
    fn test_file_not_found() {
        let result = parse_file("non_existent_file.txt");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Failed to read file"));
    }
}
