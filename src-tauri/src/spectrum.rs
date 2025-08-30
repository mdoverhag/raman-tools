use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

#[derive(Debug, Serialize, Deserialize)]
pub struct Spectrum {
    pub filename: String,
    pub filepath: String,
    pub wavenumbers: Vec<f32>,
    pub intensities: Vec<f32>,
}

pub fn parse_files(filepaths: Vec<String>) -> Result<Vec<Spectrum>, String> {
    let mut spectra = Vec::new();

    for filepath in filepaths {
        let path = Path::new(&filepath);
        let filename = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        let content = fs::read_to_string(&filepath)
            .map_err(|e| format!("Failed to read file {}: {}", filepath, e))?;

        let mut wavenumbers = Vec::new();
        let mut intensities = Vec::new();

        for line in content.lines() {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 2 {
                if let (Ok(wn), Ok(intensity)) = (parts[0].parse::<f32>(), parts[1].parse::<f32>())
                {
                    wavenumbers.push(wn);
                    intensities.push(intensity);
                }
            }
        }

        if !wavenumbers.is_empty() {
            spectra.push(Spectrum {
                filename,
                filepath,
                wavenumbers,
                intensities,
            });
        }
    }

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

        let result = parse_files(vec![file_path.to_str().unwrap().to_string()]);

        assert!(result.is_ok());
        let spectra = result.unwrap();
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].filename, "test_spectrum.txt");
        assert_eq!(spectra[0].wavenumbers, vec![200.0, 201.0, 202.0, 203.0]);
        assert_eq!(spectra[0].intensities, vec![145.0, 143.0, 143.0, 144.0]);
    }

    #[test]
    fn test_parse_multiple_files() {
        let temp_dir = TempDir::new().unwrap();
        let content1 = "200\t145\n201\t143";
        let content2 = "200\t150\n201\t148";
        let file1 = create_test_spectrum_file(&temp_dir, "spectrum1.txt", content1);
        let file2 = create_test_spectrum_file(&temp_dir, "spectrum2.txt", content2);

        let result = parse_files(vec![
            file1.to_str().unwrap().to_string(),
            file2.to_str().unwrap().to_string(),
        ]);

        assert!(result.is_ok());
        let spectra = result.unwrap();
        assert_eq!(spectra.len(), 2);
    }

    #[test]
    fn test_parse_empty_file() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = create_test_spectrum_file(&temp_dir, "empty.txt", "");

        let result = parse_files(vec![file_path.to_str().unwrap().to_string()]);

        assert!(result.is_ok());
        let spectra = result.unwrap();
        assert_eq!(spectra.len(), 0);
    }

    #[test]
    fn test_parse_malformed_lines() {
        let temp_dir = TempDir::new().unwrap();
        let content = "200\t145\ninvalid line\n202\t143\n203\tnot_a_number";
        let file_path = create_test_spectrum_file(&temp_dir, "malformed.txt", content);

        let result = parse_files(vec![file_path.to_str().unwrap().to_string()]);

        assert!(result.is_ok());
        let spectra = result.unwrap();
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].wavenumbers.len(), 2); // Only valid lines
        assert_eq!(spectra[0].intensities.len(), 2);
    }

    #[test]
    fn test_file_not_found() {
        let result = parse_files(vec!["non_existent_file.txt".to_string()]);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Failed to read file"));
    }
}
