mod python_bridge;
mod python_setup;
mod spectrum;
mod uv_installer;

use spectrum::Spectrum;

#[tauri::command]
fn parse_spectrum_files(filepaths: Vec<String>) -> Result<Vec<Spectrum>, String> {
    spectrum::parse_files(filepaths)
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
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            parse_spectrum_files,
            apply_baseline_correction,
            uv_installer::check_uv_status,
            uv_installer::download_uv,
            python_setup::check_python_status,
            python_setup::setup_python_env
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
