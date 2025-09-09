use crate::batch_baseline::{process_batch_streaming, BaselineParams, BatchUpdate};
use crate::samples::SampleStorage;
use crate::spectrum::{Spectrum, SpectrumStorage};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;
use tauri::{AppHandle, Emitter};
use uuid::Uuid;

// Helper function to convert f32 to u16 with validation
fn f32_to_u16(value: f32) -> Result<u16, String> {
    if value < 0.0 {
        return Err(format!("negative value ({})", value));
    }
    if value > 65535.0 {
        return Err(format!("value ({}) exceeds maximum (65535)", value));
    }
    if value.fract() != 0.0 {
        return Err(format!("fractional value ({})", value));
    }
    Ok(value as u16)
}

fn validate_wavenumber_sequence(wavenumbers: &[u16]) -> Result<(u16, u16, u16), String> {
    if wavenumbers.is_empty() {
        return Err("No wavenumbers found".to_string());
    }

    let first = wavenumbers[0];
    let last = *wavenumbers.last().unwrap();

    // For single point, use step of 1
    if wavenumbers.len() == 1 {
        return Ok((first, last, 1));
    }

    // Calculate expected step from first two points
    let first_step = wavenumbers[1].saturating_sub(wavenumbers[0]);
    if first_step == 0 {
        return Err("wavenumbers must be in ascending order".to_string());
    }

    // Validate that all steps are consistent
    for i in 2..wavenumbers.len() {
        let step = wavenumbers[i].saturating_sub(wavenumbers[i - 1]);
        if step != first_step {
            return Err(format!(
                "inconsistent wavenumber steps: expected {}, found {} at position {}",
                first_step, step, i
            ));
        }
    }

    Ok((first, last, first_step))
}

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

    let mut wavenumbers: Vec<u16> = Vec::new();
    let mut intensities: Vec<u16> = Vec::new();

    for (line_num, line) in content.lines().enumerate() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() == 2 {
            if let (Ok(wn_f32), Ok(intensity_f32)) =
                (parts[0].parse::<f32>(), parts[1].parse::<f32>())
            {
                // Convert wavenumber to u16
                let wn = f32_to_u16(wn_f32).map_err(|e| {
                    format!(
                        "Invalid wavenumber in {} at line {}: {}",
                        filename,
                        line_num + 1,
                        e
                    )
                })?;

                // Convert intensity to u16
                let intensity = f32_to_u16(intensity_f32).map_err(|e| {
                    format!(
                        "Invalid intensity in {} at line {}: {}",
                        filename,
                        line_num + 1,
                        e
                    )
                })?;

                wavenumbers.push(wn);
                intensities.push(intensity);
            }
        }
    }

    if wavenumbers.is_empty() {
        return Err(format!("No valid spectrum data found in {}", filename));
    }

    // Validate wavenumbers form a consistent sequence and extract parameters
    let (wavenumber_start, wavenumber_end, wavenumber_step) =
        validate_wavenumber_sequence(&wavenumbers)
            .map_err(|e| format!("Invalid wavenumber sequence in {}: {}", filename, e))?;

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
    spectrum_storage: &SpectrumStorage,
    sample_storage: &SampleStorage,
    sample_id: Uuid,
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
                // Save spectrum to storage immediately after parsing
                let saved_spectrum = spectrum_storage.create_spectrum(spectrum)?;

                // Add spectrum ID to the sample immediately
                sample_storage.add_spectra_to_sample(sample_id, vec![saved_spectrum.id])?;

                // Emit spectrum ready event with saved spectrum (includes ID)
                app.emit(
                    "import:spectrum_ready",
                    ImportEvent::SpectrumReady {
                        spectrum: saved_spectrum.clone(),
                    },
                )
                .map_err(|e| format!("Failed to emit spectrum ready event: {}", e))?;

                spectra.push(saved_spectrum);
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

                                // Update the spectrum in storage with baseline data
                                let updated_spectrum = spectrum_storage.update_spectrum(
                                    spectrum.id,
                                    spectrum_result.baseline.clone(),
                                    spectrum_result.corrected.clone(),
                                )?;

                                // Update local spectrum with the updated data
                                spectrum.baseline = Some(spectrum_result.baseline);
                                spectrum.corrected = Some(spectrum_result.corrected);

                                // Emit spectrum updated event immediately
                                app.emit(
                                    "import:spectrum_updated",
                                    ImportEvent::SpectrumUpdated {
                                        spectrum: updated_spectrum,
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

    #[test]
    fn test_negative_intensity_values() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\n201\t-50\n202\t143";
        let file_path = create_test_spectrum_file(&temp_dir, "negative.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("negative value"));
        assert!(error.contains("line 2"));
        assert!(error.contains("-50"));
    }

    #[test]
    fn test_intensity_exceeds_maximum() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\n201\t70000\n202\t143";
        let file_path = create_test_spectrum_file(&temp_dir, "overflow.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("exceeds maximum"));
        assert!(error.contains("line 2"));
        assert!(error.contains("70000"));
    }

    #[test]
    fn test_fractional_intensity_values() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\n201\t143.5\n202\t143";
        let file_path = create_test_spectrum_file(&temp_dir, "fractional.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("fractional value"));
        assert!(error.contains("line 2"));
        assert!(error.contains("143.5"));
    }

    #[test]
    fn test_valid_edge_cases() {
        let temp_dir = TempDir::new().unwrap();
        // Test edge values: 0 and 65535
        let content = "200\t0\n201\t65535\n202\t1000";
        let file_path = create_test_spectrum_file(&temp_dir, "edge_cases.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_ok());
        let spectrum = result.unwrap();
        assert_eq!(spectrum.intensities, vec![0, 65535, 1000]);
    }

    #[test]
    fn test_wavenumber_negative() {
        let temp_dir = TempDir::new().unwrap();
        let content = "-100\t145\n201\t143\n202\t143";
        let file_path = create_test_spectrum_file(&temp_dir, "negative_wn.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("negative value"));
        assert!(error.contains("-100"));
    }

    #[test]
    fn test_wavenumber_too_large() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\n201\t143\n70000\t143";
        let file_path = create_test_spectrum_file(&temp_dir, "large_wn.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("exceeds maximum"));
        assert!(error.contains("70000"));
    }

    #[test]
    fn test_descending_wavenumbers() {
        let temp_dir = TempDir::new().unwrap();
        let content = "203\t145\n202\t143\n201\t143\n200\t144";
        let file_path = create_test_spectrum_file(&temp_dir, "descending.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        // With u16 arithmetic, 202 - 203 saturates to 0
        assert!(error.contains("ascending order"));
    }

    #[test]
    fn test_zero_step_wavenumbers() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\n200\t143\n201\t143";
        let file_path = create_test_spectrum_file(&temp_dir, "zero_step.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("ascending order"));
    }

    #[test]
    fn test_fractional_wavenumber_rejected() {
        let temp_dir = TempDir::new().unwrap();
        // Fractional wavenumbers should be rejected
        let content = "200.4\t145\n201.4\t143\n202.4\t143";
        let file_path = create_test_spectrum_file(&temp_dir, "fractional_wn.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("fractional value"));
        assert!(error.contains("200.4"));
    }

    #[test]
    fn test_wavenumber_at_boundary() {
        let temp_dir = TempDir::new().unwrap();
        // Test wavenumbers at the exact boundaries - need consistent step!
        let content = "0\t100\n1\t101\n2\t102";
        let file_path = create_test_spectrum_file(&temp_dir, "boundary.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_ok());
        let spectrum = result.unwrap();
        assert_eq!(spectrum.wavenumber_start, 0);
        assert_eq!(spectrum.wavenumber_end, 2);
        assert_eq!(spectrum.wavenumber_step, 1);

        // Test upper boundary
        let content = "65533\t100\n65534\t101\n65535\t102";
        let file_path = create_test_spectrum_file(&temp_dir, "boundary_upper.txt", content);

        let result = parse_file(file_path.to_str().unwrap());
        assert!(result.is_ok());
        let spectrum = result.unwrap();
        assert_eq!(spectrum.wavenumber_start, 65533);
        assert_eq!(spectrum.wavenumber_end, 65535);
    }

    #[test]
    fn test_inconsistent_wavenumber_steps() {
        let temp_dir = TempDir::new().unwrap();
        // Wavenumbers with inconsistent steps: 1, 1, 2
        let content = "200\t145\n201\t143\n202\t144\n204\t146";
        let file_path = create_test_spectrum_file(&temp_dir, "inconsistent.txt", content);

        let result = parse_file(file_path.to_str().unwrap());

        assert!(result.is_err());
        let error = result.unwrap_err();
        assert!(error.contains("inconsistent wavenumber steps"));
        assert!(error.contains("position 3"));
    }
}
