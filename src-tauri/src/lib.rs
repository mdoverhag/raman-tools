mod python_bridge;
mod python_setup;
mod samples;
mod spectrum;
mod uv_installer;

use samples::{Sample, SampleStorage, UpdateSampleData};
use spectrum::Spectrum;
use tauri::State;
use uuid::Uuid;

#[tauri::command]
async fn parse_spectrum_files(
    app: tauri::AppHandle,
    filepaths: Vec<String>,
) -> Result<Vec<Spectrum>, String> {
    // Parse the files first
    let mut spectra = spectrum::parse_files(filepaths)?;

    // Apply baseline correction to each spectrum
    // Using default parameters for now - could be made configurable later
    let denoise = true;
    let window_size = 5;
    let lambda_param = 1e7;
    let p = 0.01;
    let d = 2;

    for spectrum in &mut spectra {
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
            }
            Err(e) => {
                eprintln!(
                    "Failed to apply baseline correction to {}: {}",
                    spectrum.filename, e
                );
                // Continue without baseline for this spectrum
            }
        }
    }

    Ok(spectra)
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
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            parse_spectrum_files,
            create_sample,
            list_samples,
            get_sample,
            update_sample,
            delete_sample,
            apply_baseline_correction,
            uv_installer::check_uv_status,
            uv_installer::download_uv,
            python_setup::check_python_status,
            python_setup::setup_python_env
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
