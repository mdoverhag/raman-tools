mod python_bridge;
mod spectrum;

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

#[tauri::command]
fn check_python_runtime(app: tauri::AppHandle) -> Result<bool, String> {
    python_bridge::check_python_availability(&app)
}

#[tauri::command]
fn get_python_info(app: tauri::AppHandle) -> Result<String, String> {
    python_bridge::get_python_info(&app)
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            parse_spectrum_files,
            apply_baseline_correction,
            check_python_runtime,
            get_python_info
        ])
        .setup(|app| {
            // Check Python runtime availability on startup (development only)
            #[cfg(debug_assertions)]
            {
                let app_handle = app.handle();

                match python_bridge::check_python_availability(&app_handle) {
                    Ok(true) => {
                        println!("✅ Python runtime found");

                        // Test baseline correction with a simple spectrum
                        println!("   Testing baseline correction...");
                        let test_spectrum = vec![1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];

                        // Use Tauri's built-in async runtime to block on the async function
                        let app_handle_clone = app_handle.clone();

                        match tauri::async_runtime::block_on(python_bridge::apply_baseline_correction(
                            app_handle_clone,
                            test_spectrum,
                            false,  // denoise
                            5,      // window_size
                            1e7,    // lambda_param
                            0.01,   // p
                            2,      // d
                        )) {
                            Ok(_) => {
                                println!("✅ Baseline correction test successful - all systems ready!");
                            }
                            Err(e) => {
                                eprintln!("\n❌ Baseline correction test failed: {}", e);
                                eprintln!("\nThe Python runtime was found but baseline correction is not working properly.");
                                eprintln!("This may indicate missing dependencies (numpy/scipy).");
                                eprintln!("\nPlease rebuild the Python bundle:\n");
                                eprintln!("    cd src-tauri && ./build-python.sh\n");

                                // Exit cleanly without panic, this should only happen in development
                                std::process::exit(1);
                            }
                        }
                    }
                    Ok(false) | Err(_) => {
                        eprintln!("\n❌ Python runtime not found - cannot start application");
                        eprintln!("\nBaseline correction requires a bundled Python runtime with numpy and scipy.");
                        eprintln!("Please run the following command to set up the Python bundle:\n");
                        eprintln!("    cd src-tauri && ./build-python.sh\n");
                        eprintln!("Then restart the application.\n");

                        // Exit cleanly without panic, this should only happen in development
                        std::process::exit(1);
                    }
                }
            }

            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
