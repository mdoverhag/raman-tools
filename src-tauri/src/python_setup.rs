//! Python Environment Setup
//!
//! Uses the installed uv to create Python virtual environment and install packages

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tauri::{AppHandle, Emitter, Manager};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SetupProgress {
    pub stage: String,
    pub message: String,
    pub percentage: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PythonStatus {
    pub installed: bool,
    pub python_version: Option<String>,
    pub numpy_version: Option<String>,
    pub scipy_version: Option<String>,
}

/// Get the runtime directory
pub fn get_runtime_dir(app: &AppHandle) -> Result<PathBuf, String> {
    let data_dir = app
        .path()
        .app_data_dir()
        .map_err(|e| format!("Failed to get app data directory: {}", e))?;

    Ok(data_dir.join("runtime"))
}

/// Get the path to uv
fn get_uv_path(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_dir = get_runtime_dir(app)?;
    let uv_name = if cfg!(target_os = "windows") {
        "uv.exe"
    } else {
        "uv"
    };

    let path = runtime_dir.join("bin").join(uv_name);
    if !path.exists() {
        return Err("uv is not installed".to_string());
    }

    Ok(path)
}

/// Get the path to the virtual environment
fn get_venv_dir(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_dir = get_runtime_dir(app)?;
    Ok(runtime_dir.join("venv"))
}

/// Get the path to Python in the venv
pub fn get_python_path(app: &AppHandle) -> Result<PathBuf, String> {
    let venv_dir = get_venv_dir(app)?;

    let python_path = if cfg!(target_os = "windows") {
        venv_dir.join("Scripts").join("python.exe")
    } else {
        venv_dir.join("bin").join("python")
    };

    if !python_path.exists() {
        return Err("Python is not installed".to_string());
    }

    Ok(python_path)
}

/// Get the path to the baseline correction script
pub fn get_baseline_script_path(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_dir = get_runtime_dir(app)?;
    let script_path = runtime_dir.join("baseline_correction.py");

    if !script_path.exists() {
        return Err("Baseline correction script not found".to_string());
    }

    Ok(script_path)
}

/// Get the path to the batch processor script
pub fn get_batch_processor_path(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_dir = get_runtime_dir(app)?;
    let script_path = runtime_dir.join("batch_processor.py");

    if !script_path.exists() {
        return Err("Batch processor script not found".to_string());
    }

    Ok(script_path)
}

/// Get the path to the calc averages script
pub fn get_calc_averages_path(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_dir = get_runtime_dir(app)?;
    let script_path = runtime_dir.join("calc_averages.py");

    if !script_path.exists() {
        return Err("Calculate averages script not found".to_string());
    }

    Ok(script_path)
}

/// Sync Python files to runtime directory
/// This ensures all embedded Python scripts are copied to the runtime directory
/// Called on every app start to handle updates
pub fn sync_python_files(app: AppHandle) -> Result<(), String> {
    let runtime_dir = get_runtime_dir(&app)?;

    // Create runtime directory if it doesn't exist
    if !runtime_dir.exists() {
        fs::create_dir_all(&runtime_dir)
            .map_err(|e| format!("Failed to create runtime directory: {}", e))?;
    }

    // Copy/update baseline correction script
    let baseline_content = include_str!("../python/baseline_correction.py");
    let baseline_path = runtime_dir.join("baseline_correction.py");
    fs::write(&baseline_path, baseline_content)
        .map_err(|e| format!("Failed to write baseline script: {}", e))?;

    // Copy/update batch processor script
    let batch_content = include_str!("../python/batch_processor.py");
    let batch_path = runtime_dir.join("batch_processor.py");
    fs::write(&batch_path, batch_content)
        .map_err(|e| format!("Failed to write batch processor script: {}", e))?;

    // Copy/update calc averages script
    let calc_averages_content = include_str!("../python/calc_averages.py");
    let calc_averages_path = runtime_dir.join("calc_averages.py");
    fs::write(&calc_averages_path, calc_averages_content)
        .map_err(|e| format!("Failed to write calc averages script: {}", e))?;

    // Copy/update normalize spectra script
    let normalize_content = include_str!("../python/normalize_spectra.py");
    let normalize_path = runtime_dir.join("normalize_spectra.py");
    fs::write(&normalize_path, normalize_content)
        .map_err(|e| format!("Failed to write normalize spectra script: {}", e))?;

    // Copy/update deconvolute NNLS script
    let deconvolute_content = include_str!("../python/deconvolute_nnls.py");
    let deconvolute_path = runtime_dir.join("deconvolute_nnls.py");
    fs::write(&deconvolute_path, deconvolute_content)
        .map_err(|e| format!("Failed to write deconvolute NNLS script: {}", e))?;

    // Copy/update requirements.txt
    let requirements_content = include_str!("../python/requirements.txt");
    let requirements_path = runtime_dir.join("requirements.txt");
    fs::write(&requirements_path, requirements_content)
        .map_err(|e| format!("Failed to write requirements: {}", e))?;

    Ok(())
}

/// Check if Python environment is set up
#[tauri::command]
pub fn check_python_status(app: AppHandle) -> PythonStatus {
    let python_path = match get_python_path(&app) {
        Ok(path) => path,
        Err(_) => {
            return PythonStatus {
                installed: false,
                python_version: None,
                numpy_version: None,
                scipy_version: None,
            }
        }
    };

    // Get Python version
    let python_version = Command::new(&python_path)
        .arg("--version")
        .output()
        .ok()
        .and_then(|output| {
            String::from_utf8(output.stdout)
                .ok()
                .map(|s| s.trim().replace("Python ", ""))
        });

    // Check numpy
    let numpy_version = Command::new(&python_path)
        .args(["-c", "import numpy; print(numpy.__version__)"])
        .output()
        .ok()
        .and_then(|output| {
            if output.status.success() {
                String::from_utf8(output.stdout)
                    .ok()
                    .map(|s| s.trim().to_string())
            } else {
                None
            }
        });

    // Check scipy
    let scipy_version = Command::new(&python_path)
        .args(["-c", "import scipy; print(scipy.__version__)"])
        .output()
        .ok()
        .and_then(|output| {
            if output.status.success() {
                String::from_utf8(output.stdout)
                    .ok()
                    .map(|s| s.trim().to_string())
            } else {
                None
            }
        });

    PythonStatus {
        installed: python_version.is_some(),
        python_version,
        numpy_version,
        scipy_version,
    }
}

/// Set up Python environment
#[tauri::command]
pub async fn setup_python_env(app: AppHandle) -> Result<(), String> {
    let uv_path = get_uv_path(&app)?;
    let venv_dir = get_venv_dir(&app)?;
    let runtime_dir = get_runtime_dir(&app)?;

    // Emit progress: Creating venv
    app.emit(
        "python-setup-progress",
        SetupProgress {
            stage: "creating_venv".to_string(),
            message: "Creating Python virtual environment...".to_string(),
            percentage: Some(20),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Create virtual environment
    let output = Command::new(&uv_path)
        .args(["venv", "--python", "3.13"])
        .arg(&venv_dir)
        .output()
        .map_err(|e| format!("Failed to create venv: {}", e))?;

    if !output.status.success() {
        return Err(format!(
            "Failed to create venv: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }

    // Requirements.txt is already synced via sync_python_files()
    let requirements_path = runtime_dir.join("requirements.txt");

    // Emit progress: Installing packages
    app.emit(
        "python-setup-progress",
        SetupProgress {
            stage: "installing_packages".to_string(),
            message: "Installing scientific packages (numpy, scipy)...".to_string(),
            percentage: Some(40),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Install from requirements.txt
    let output = Command::new(&uv_path)
        .args(["pip", "install", "--python"])
        .arg(&venv_dir)
        .arg("-r")
        .arg(&requirements_path)
        .output()
        .map_err(|e| format!("Failed to install packages: {}", e))?;

    if !output.status.success() {
        return Err(format!(
            "Failed to install packages: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }

    // Emit progress: Copying script
    app.emit(
        "python-setup-progress",
        SetupProgress {
            stage: "copying_script".to_string(),
            message: "Installing baseline correction script...".to_string(),
            percentage: Some(70),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Python files are now synced on app start via sync_python_files()
    // No need to copy them here during initial setup

    // Emit progress: Testing
    app.emit(
        "python-setup-progress",
        SetupProgress {
            stage: "testing".to_string(),
            message: "Testing Python environment...".to_string(),
            percentage: Some(85),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Test that everything works
    let python_path = get_python_path(&app)?;
    let test_script = r#"
import numpy as np
import scipy
print("OK")
"#;

    let output = Command::new(&python_path)
        .args(["-c", test_script])
        .output()
        .map_err(|e| format!("Failed to test Python: {}", e))?;

    if !output.status.success() {
        return Err(format!(
            "Python test failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }

    // Emit progress: Complete
    app.emit(
        "python-setup-progress",
        SetupProgress {
            stage: "complete".to_string(),
            message: "Python environment ready!".to_string(),
            percentage: Some(100),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    Ok(())
}
