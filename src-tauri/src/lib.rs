mod batch_baseline;
mod python_bridge;
mod python_setup;
mod samples;
mod spectrum;
mod spectrum_importer;
mod uv_installer;

use samples::{Sample, SampleStorage, UpdateSampleData};
use spectrum::{Spectrum, SpectrumStorage};
use tauri::State;
use uuid::Uuid;

#[tauri::command]
async fn parse_spectrum_files(
    app: tauri::AppHandle,
    filepaths: Vec<String>,
    sample_id: String, // Now required since drop zone only shows with selected sample
    sample_storage: State<'_, SampleStorage>,
    spectrum_storage: State<'_, SpectrumStorage>,
) -> Result<Vec<Spectrum>, String> {
    // Parse the sample ID
    let sample_uuid =
        Uuid::parse_str(&sample_id).map_err(|e| format!("Invalid sample ID: {}", e))?;

    // Import the spectra
    let spectra = spectrum_importer::import_spectra(app, filepaths).await?;

    // Save each spectrum to storage
    let mut saved_spectra = Vec::new();
    for spectrum in spectra {
        let saved = spectrum_storage.create_spectrum(spectrum)?;
        saved_spectra.push(saved);
    }

    // Collect the spectrum IDs
    let spectrum_ids: Vec<Uuid> = saved_spectra.iter().map(|s| s.id).collect();

    // Add spectrum IDs to the sample
    sample_storage.add_spectra_to_sample(sample_uuid, spectrum_ids.clone())?;

    Ok(saved_spectra)
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
fn delete_sample(storage: State<SampleStorage>, id: Uuid) -> Result<(), String> {
    storage.delete_sample(id)
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

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .manage(SampleStorage::new())
        .manage(SpectrumStorage::new())
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
