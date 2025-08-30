mod spectrum;

use spectrum::Spectrum;

#[tauri::command]
fn parse_spectrum_files(filepaths: Vec<String>) -> Result<Vec<Spectrum>, String> {
    spectrum::parse_files(filepaths)
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![parse_spectrum_files])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
