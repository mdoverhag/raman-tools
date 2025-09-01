//! UV Package Manager Installer
//!
//! Downloads and manages the uv Python package manager in app data directory

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use tauri::{AppHandle, Emitter, Manager};

const UV_VERSION: &str = "0.8.14";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UvStatus {
    pub installed: bool,
    pub version: Option<String>,
    pub path: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DownloadProgress {
    pub stage: String,
    pub message: String,
    pub percentage: Option<u32>,
}

/// Get the runtime directory in app data
fn get_runtime_dir(app: &AppHandle) -> Result<PathBuf, String> {
    let data_dir = app
        .path()
        .app_data_dir()
        .map_err(|e| format!("Failed to get app data directory: {}", e))?;

    Ok(data_dir.join("runtime"))
}

/// Get the path where uv should be installed
fn get_uv_path(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_dir = get_runtime_dir(app)?;
    let uv_name = if cfg!(target_os = "windows") {
        "uv.exe"
    } else {
        "uv"
    };

    Ok(runtime_dir.join("bin").join(uv_name))
}

/// Check if uv is installed and get its version
#[tauri::command]
pub fn check_uv_status(app: AppHandle) -> UvStatus {
    let uv_path = match get_uv_path(&app) {
        Ok(path) => path,
        Err(_) => {
            return UvStatus {
                installed: false,
                version: None,
                path: None,
            }
        }
    };

    if !uv_path.exists() {
        return UvStatus {
            installed: false,
            version: None,
            path: None,
        };
    }

    // Try to get version
    let version = Command::new(&uv_path)
        .arg("--version")
        .output()
        .ok()
        .and_then(|output| {
            if output.status.success() {
                String::from_utf8(output.stdout)
                    .ok()
                    .map(|s| s.trim().replace("uv ", ""))
            } else {
                None
            }
        });

    UvStatus {
        installed: true,
        version,
        path: Some(uv_path.to_string_lossy().to_string()),
    }
}

/// Download and install uv
#[tauri::command]
pub async fn download_uv(app: AppHandle) -> Result<(), String> {
    let runtime_dir = get_runtime_dir(&app)?;
    let bin_dir = runtime_dir.join("bin");

    // Create directories
    fs::create_dir_all(&bin_dir).map_err(|e| format!("Failed to create directories: {}", e))?;

    // Emit progress: Starting
    app.emit(
        "uv-download-progress",
        DownloadProgress {
            stage: "starting".to_string(),
            message: "Preparing to download uv...".to_string(),
            percentage: Some(0),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Determine platform and download URL
    let download_url = if cfg!(target_os = "windows") {
        format!(
            "https://github.com/astral-sh/uv/releases/download/{}/uv-x86_64-pc-windows-msvc.zip",
            UV_VERSION
        )
    } else if cfg!(target_os = "macos") {
        // Check architecture
        let arch = if cfg!(target_arch = "aarch64") {
            "aarch64"
        } else {
            "x86_64"
        };
        format!(
            "https://github.com/astral-sh/uv/releases/download/{}/uv-{}-apple-darwin.tar.gz",
            UV_VERSION, arch
        )
    } else if cfg!(target_os = "linux") {
        format!("https://github.com/astral-sh/uv/releases/download/{}/uv-x86_64-unknown-linux-gnu.tar.gz", UV_VERSION)
    } else {
        return Err("Unsupported platform".to_string());
    };

    // Emit progress: Downloading
    app.emit(
        "uv-download-progress",
        DownloadProgress {
            stage: "downloading".to_string(),
            message: format!("Downloading uv {}...", UV_VERSION),
            percentage: Some(30),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Download the file
    let response = reqwest::get(&download_url)
        .await
        .map_err(|e| format!("Failed to download uv: {}", e))?;

    if !response.status().is_success() {
        return Err(format!(
            "Download failed with status: {}",
            response.status()
        ));
    }

    let bytes = response
        .bytes()
        .await
        .map_err(|e| format!("Failed to read download: {}", e))?;

    // Emit progress: Extracting
    app.emit(
        "uv-download-progress",
        DownloadProgress {
            stage: "extracting".to_string(),
            message: "Extracting uv...".to_string(),
            percentage: Some(60),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Save and extract based on platform
    if cfg!(target_os = "windows") {
        // Save zip file
        let zip_path = runtime_dir.join("uv.zip");
        fs::write(&zip_path, &bytes).map_err(|e| format!("Failed to write zip file: {}", e))?;

        // Extract using PowerShell
        let output = Command::new("powershell")
            .args([
                "-Command",
                &format!(
                    "Expand-Archive -Path '{}' -DestinationPath '{}' -Force; Get-ChildItem '{}' -Recurse -Filter 'uv.exe' | Select-Object -First 1 | Move-Item -Destination '{}' -Force",
                    zip_path.display(),
                    runtime_dir.display(),
                    runtime_dir.display(),
                    bin_dir.join("uv.exe").display()
                )
            ])
            .output()
            .map_err(|e| format!("Failed to extract uv: {}", e))?;

        if !output.status.success() {
            return Err(format!(
                "Failed to extract uv: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }

        // Clean up
        let _ = fs::remove_file(zip_path);
    } else {
        // Save tar.gz file
        let tar_path = runtime_dir.join("uv.tar.gz");
        fs::write(&tar_path, &bytes).map_err(|e| format!("Failed to write tar file: {}", e))?;

        // Extract to temp directory first
        let temp_dir = runtime_dir.join("temp");
        fs::create_dir_all(&temp_dir).map_err(|e| format!("Failed to create temp dir: {}", e))?;

        // Extract using tar
        let output = Command::new("tar")
            .args([
                "-xzf",
                tar_path.to_str().unwrap(),
                "-C",
                temp_dir.to_str().unwrap(),
            ])
            .output()
            .map_err(|e| format!("Failed to extract uv: {}", e))?;

        if !output.status.success() {
            return Err(format!(
                "Failed to extract uv: {}",
                String::from_utf8_lossy(&output.stderr)
            ));
        }

        // Find the extracted uv binary and move it to bin
        // The archive contains a directory like "uv-aarch64-apple-darwin/" with uv and uvx inside
        let entries =
            fs::read_dir(&temp_dir).map_err(|e| format!("Failed to read temp dir: {}", e))?;

        for entry in entries {
            let entry = entry.map_err(|e| format!("Failed to read dir entry: {}", e))?;
            let path = entry.path();

            if path.is_dir() {
                // This should be the uv-* directory
                let uv_src = path.join("uv");
                if uv_src.exists() {
                    let uv_dst = bin_dir.join("uv");
                    fs::rename(&uv_src, &uv_dst)
                        .map_err(|e| format!("Failed to move uv binary: {}", e))?;

                    // Also move uvx if it exists
                    let uvx_src = path.join("uvx");
                    if uvx_src.exists() {
                        let uvx_dst = bin_dir.join("uvx");
                        let _ = fs::rename(&uvx_src, &uvx_dst);
                    }

                    break;
                }
            }
        }

        // Make executable
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let uv_path = bin_dir.join("uv");
            let mut perms = fs::metadata(&uv_path)
                .map_err(|e| format!("Failed to get file metadata: {}", e))?
                .permissions();
            perms.set_mode(0o755);
            fs::set_permissions(&uv_path, perms)
                .map_err(|e| format!("Failed to set permissions: {}", e))?;

            // Also make uvx executable if it exists
            let uvx_path = bin_dir.join("uvx");
            if uvx_path.exists() {
                let mut perms = fs::metadata(&uvx_path)
                    .map_err(|e| format!("Failed to get uvx metadata: {}", e))?
                    .permissions();
                perms.set_mode(0o755);
                let _ = fs::set_permissions(&uvx_path, perms);
            }
        }

        // Clean up
        let _ = fs::remove_file(tar_path);
        let _ = fs::remove_dir_all(temp_dir);
    }

    // Emit progress: Verifying
    app.emit(
        "uv-download-progress",
        DownloadProgress {
            stage: "verifying".to_string(),
            message: "Verifying installation...".to_string(),
            percentage: Some(90),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    // Verify installation
    let uv_path = get_uv_path(&app)?;
    let output = Command::new(&uv_path)
        .arg("--version")
        .output()
        .map_err(|e| format!("Failed to run uv: {}", e))?;

    if !output.status.success() {
        return Err("Failed to verify uv installation".to_string());
    }

    let version = String::from_utf8_lossy(&output.stdout);

    // Emit progress: Complete
    app.emit(
        "uv-download-progress",
        DownloadProgress {
            stage: "complete".to_string(),
            message: format!("Successfully installed {}", version.trim()),
            percentage: Some(100),
        },
    )
    .map_err(|e| format!("Failed to emit progress: {}", e))?;

    Ok(())
}
