mod batch_baseline;
mod deconvolution;
mod molecules;
mod python_bridge;
mod python_setup;
mod samples;
mod spectrum;
mod spectrum_importer;
mod uv_installer;

use deconvolution::{CreateDeconvolutionData, DeconvolutionRun, DeconvolutionStorage};
use samples::{Sample, SampleStorage, UpdateSampleData};
use serde::Serialize;
use spectrum::{Spectrum, SpectrumStorage};
use std::collections::HashMap;
use tauri::State;
use uuid::Uuid;

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ImportResult {
    success: bool,
    count: usize,
    sample_id: String,
}

#[tauri::command]
async fn parse_spectrum_files(
    app: tauri::AppHandle,
    filepaths: Vec<String>,
    sample_id: String, // Now required since drop zone only shows with selected sample
    sample_storage: State<'_, SampleStorage>,
    spectrum_storage: State<'_, SpectrumStorage>,
) -> Result<ImportResult, String> {
    // Parse the sample ID
    let sample_uuid =
        Uuid::parse_str(&sample_id).map_err(|e| format!("Invalid sample ID: {}", e))?;

    // Import the spectra (they're saved to storage and linked to sample during import)
    let spectra = spectrum_importer::import_spectra(
        app,
        filepaths,
        &spectrum_storage,
        &sample_storage,
        sample_uuid,
    )
    .await?;

    let count = spectra.len();

    Ok(ImportResult {
        success: count > 0,
        count,
        sample_id,
    })
}

#[tauri::command]
fn create_sample(storage: State<SampleStorage>, name: String) -> Result<Sample, String> {
    storage.create_sample(name)
}

#[tauri::command]
fn list_samples(storage: State<SampleStorage>) -> Result<Vec<Sample>, String> {
    storage.list_samples()
}

#[tauri::command]
fn get_sample(storage: State<SampleStorage>, id: Uuid) -> Result<Sample, String> {
    storage.get_sample(id)
}

#[tauri::command]
fn update_sample(
    storage: State<SampleStorage>,
    id: Uuid,
    updates: UpdateSampleData,
) -> Result<Sample, String> {
    storage.update_sample(id, updates)
}

#[tauri::command]
fn delete_sample(
    sample_storage: State<SampleStorage>,
    spectrum_storage: State<SpectrumStorage>,
    id: Uuid,
) -> Result<(), String> {
    // First get the sample to find its spectrum IDs
    let sample = sample_storage.get_sample(id)?;

    // Delete all spectra associated with this sample
    for spectrum_id in sample.spectrum_ids {
        // Ignore errors for missing spectra (they may have been deleted already)
        let _ = spectrum_storage.delete_spectrum(spectrum_id);
    }

    // Now delete the sample
    sample_storage.delete_sample(id)
}

#[tauri::command]
fn get_spectrum(storage: State<SpectrumStorage>, id: Uuid) -> Result<Spectrum, String> {
    storage.get_spectrum(id)
}

#[tauri::command]
async fn apply_baseline_correction(
    app: tauri::AppHandle,
    spectrum: Vec<f64>,
    denoise: bool,
    window_size: i32,
    lambda_param: f64,
    p: f64,
    d: i32,
) -> Result<python_bridge::BaselineResponse, String> {
    python_bridge::apply_baseline_correction(
        app,
        spectrum,
        denoise,
        window_size,
        lambda_param,
        p,
        d,
    )
    .await
}

// Deconvolution commands
#[tauri::command]
fn create_deconvolution_run(
    data: CreateDeconvolutionData,
    deconvolution_storage: State<'_, DeconvolutionStorage>,
) -> Result<DeconvolutionRun, String> {
    deconvolution_storage.create_run(
        data.multiplex_sample_id,
        data.reference_sample_ids,
        data.wavenumber_range,
    )
}

#[tauri::command]
fn get_deconvolution_run(
    id: String,
    deconvolution_storage: State<'_, DeconvolutionStorage>,
) -> Result<DeconvolutionRun, String> {
    let uuid = Uuid::parse_str(&id).map_err(|e| e.to_string())?;
    deconvolution_storage.get_run(&uuid)
}

#[tauri::command]
async fn calculate_normalization(
    app: tauri::AppHandle,
    deconvolution_id: String,
    sample_storage: State<'_, SampleStorage>,
    _spectrum_storage: State<'_, SpectrumStorage>,
    deconvolution_storage: State<'_, DeconvolutionStorage>,
) -> Result<DeconvolutionRun, String> {
    let uuid = Uuid::parse_str(&deconvolution_id).map_err(|e| e.to_string())?;
    let run = deconvolution_storage.get_run(&uuid)?;

    // Get multiplex sample and its average spectrum
    let multiplex_sample = sample_storage.get_sample(run.multiplex_sample_id)?;
    let multiplex_spectrum = multiplex_sample
        .average_corrected
        .ok_or("Multiplex sample has no average spectrum")?;

    // Get reference samples and their average spectra
    let mut reference_spectra = HashMap::new();
    for ref_id in &run.reference_sample_ids {
        let ref_sample = sample_storage.get_sample(*ref_id)?;
        let ref_spectrum = ref_sample.average_corrected.ok_or(format!(
            "Reference sample {} has no average spectrum",
            ref_sample.name
        ))?;

        // Get the Raman molecule name for this reference
        if let Some(molecule_pair) = ref_sample.molecule_pairs.first() {
            let raman_name = format!("{:?}", molecule_pair.raman);
            reference_spectra.insert(raman_name, ref_spectrum);
        }
    }

    // Call Python to normalize the spectra
    let (normalized_multiplex, normalized_references) = python_bridge::normalize_spectra(
        app,
        multiplex_spectrum,
        reference_spectra,
        run.wavenumber_range,
    )
    .await?;

    // Update the run with normalized data
    deconvolution_storage.update_normalized_data(
        &uuid,
        normalized_multiplex,
        normalized_references,
    )?;
    deconvolution_storage.get_run(&uuid)
}

#[tauri::command]
async fn perform_deconvolution(
    app: tauri::AppHandle,
    deconvolution_id: String,
    deconvolution_storage: State<'_, DeconvolutionStorage>,
) -> Result<python_bridge::DeconvolutionResponse, String> {
    let uuid = Uuid::parse_str(&deconvolution_id).map_err(|e| e.to_string())?;
    let run = deconvolution_storage.get_run(&uuid)?;

    // Check that we have normalized data
    let multiplex = run
        .normalized_multiplex
        .ok_or("No normalized multiplex data available")?;
    let references = run
        .normalized_references
        .ok_or("No normalized reference data available")?;

    // Perform NNLS deconvolution
    let result =
        python_bridge::perform_nnls_deconvolution(app, multiplex, references, run.wavenumber_range)
            .await?;

    Ok(result)
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .manage(SampleStorage::new())
        .manage(SpectrumStorage::new())
        .manage(DeconvolutionStorage::new())
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            parse_spectrum_files,
            create_sample,
            list_samples,
            get_sample,
            update_sample,
            delete_sample,
            get_spectrum,
            apply_baseline_correction,
            create_deconvolution_run,
            get_deconvolution_run,
            calculate_normalization,
            perform_deconvolution,
            uv_installer::check_uv_status,
            uv_installer::download_uv,
            python_setup::check_python_status,
            python_setup::setup_python_env
        ])
        .setup(|app| {
            // Sync Python files on every app start
            python_setup::sync_python_files(app.handle().clone())?;
            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
