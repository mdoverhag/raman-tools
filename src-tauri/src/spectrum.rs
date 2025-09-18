use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Mutex;
use uuid::Uuid;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Spectrum {
    pub id: Uuid,
    pub filename: String,
    pub wavenumber_start: u16,
    pub wavenumber_end: u16,
    pub wavenumber_step: u16,
    pub intensities: Vec<u16>,

    // Baseline calculated by ALS algorithm - stored to preserve analysis parameters
    pub baseline: Option<Vec<f64>>,

    // Baseline-corrected intensities (intensity - baseline)
    pub corrected: Option<Vec<f64>>,
}

impl Spectrum {
    pub fn new(
        filename: String,
        wavenumber_start: u16,
        wavenumber_end: u16,
        wavenumber_step: u16,
        intensities: Vec<u16>,
    ) -> Self {
        Self {
            id: Uuid::new_v4(),
            filename,
            wavenumber_start,
            wavenumber_end,
            wavenumber_step,
            intensities,
            baseline: None,
            corrected: None,
        }
    }
}

// Storage for spectra - similar to SampleStorage
pub struct SpectrumStorage {
    spectra: Mutex<HashMap<Uuid, Spectrum>>,
}

impl SpectrumStorage {
    pub fn new() -> Self {
        Self {
            spectra: Mutex::new(HashMap::new()),
        }
    }

    // Create a new spectrum and store it
    pub fn create_spectrum(&self, spectrum: Spectrum) -> Result<Spectrum, String> {
        let mut spectra = self.spectra.lock().map_err(|e| e.to_string())?;
        let id = spectrum.id;
        spectra.insert(id, spectrum.clone());
        Ok(spectrum)
    }

    // Get a spectrum by ID
    pub fn get_spectrum(&self, id: Uuid) -> Result<Spectrum, String> {
        let spectra = self.spectra.lock().map_err(|e| e.to_string())?;
        spectra
            .get(&id)
            .cloned()
            .ok_or_else(|| format!("Spectrum with id {} not found", id))
    }

    // Get multiple spectra by IDs
    pub fn get_spectra_by_ids(&self, ids: &[Uuid]) -> Result<Vec<Spectrum>, String> {
        let spectra = self.spectra.lock().map_err(|e| e.to_string())?;
        let mut result = Vec::new();

        for id in ids {
            if let Some(spectrum) = spectra.get(id) {
                result.push(spectrum.clone());
            }
        }

        Ok(result)
    }

    // Update a spectrum (e.g., after baseline correction)
    pub fn update_spectrum(
        &self,
        id: Uuid,
        baseline: Vec<f64>,
        corrected: Vec<f64>,
    ) -> Result<Spectrum, String> {
        let mut spectra = self.spectra.lock().map_err(|e| e.to_string())?;
        let spectrum = spectra
            .get_mut(&id)
            .ok_or_else(|| format!("Spectrum with id {} not found", id))?;

        spectrum.baseline = Some(baseline);
        spectrum.corrected = Some(corrected);

        Ok(spectrum.clone())
    }

    // Delete a spectrum by ID
    pub fn delete_spectrum(&self, id: Uuid) -> Result<(), String> {
        let mut spectra = self.spectra.lock().map_err(|e| e.to_string())?;
        spectra
            .remove(&id)
            .ok_or_else(|| format!("Spectrum with id {} not found", id))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_spectrum() -> Spectrum {
        Spectrum::new(
            "test.txt".to_string(),
            200,
            205,
            1,
            vec![100, 101, 102, 103, 104, 105],
        )
    }

    #[test]
    fn test_create_spectrum() {
        let storage = SpectrumStorage::new();
        let spectrum = create_test_spectrum();
        let id = spectrum.id;

        let result = storage.create_spectrum(spectrum);
        assert!(result.is_ok());

        let created = result.unwrap();
        assert_eq!(created.id, id);
        assert_eq!(created.filename, "test.txt");
        assert_eq!(created.intensities.len(), 6);
    }

    #[test]
    fn test_get_spectrum() {
        let storage = SpectrumStorage::new();
        let spectrum = create_test_spectrum();
        let id = spectrum.id;

        storage.create_spectrum(spectrum).unwrap();

        let result = storage.get_spectrum(id);
        assert!(result.is_ok());

        let retrieved = result.unwrap();
        assert_eq!(retrieved.id, id);
        assert_eq!(retrieved.filename, "test.txt");
    }

    #[test]
    fn test_get_nonexistent_spectrum() {
        let storage = SpectrumStorage::new();
        let id = Uuid::new_v4();

        let result = storage.get_spectrum(id);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_get_spectra_by_ids() {
        let storage = SpectrumStorage::new();

        // Create three test spectra
        let spectrum1 = create_test_spectrum();
        let spectrum2 = Spectrum::new(
            "test2.txt".to_string(),
            300,
            305,
            1,
            vec![200, 201, 202, 203, 204, 205],
        );
        let spectrum3 = Spectrum::new(
            "test3.txt".to_string(),
            400,
            405,
            1,
            vec![300, 301, 302, 303, 304, 305],
        );

        let id1 = spectrum1.id;
        let _id2 = spectrum2.id;
        let id3 = spectrum3.id;

        storage.create_spectrum(spectrum1).unwrap();
        storage.create_spectrum(spectrum2).unwrap();
        storage.create_spectrum(spectrum3).unwrap();

        // Get multiple spectra
        let result = storage.get_spectra_by_ids(&[id1, id3]);
        assert!(result.is_ok());

        let spectra = result.unwrap();
        assert_eq!(spectra.len(), 2);
        assert_eq!(spectra[0].filename, "test.txt");
        assert_eq!(spectra[1].filename, "test3.txt");
    }

    #[test]
    fn test_get_spectra_by_ids_with_missing() {
        let storage = SpectrumStorage::new();

        let spectrum = create_test_spectrum();
        let existing_id = spectrum.id;
        let missing_id = Uuid::new_v4();

        storage.create_spectrum(spectrum).unwrap();

        // Request both existing and non-existing IDs
        let result = storage.get_spectra_by_ids(&[existing_id, missing_id]);
        assert!(result.is_ok());

        let spectra = result.unwrap();
        assert_eq!(spectra.len(), 1); // Only returns found spectra
        assert_eq!(spectra[0].filename, "test.txt");
    }

    #[test]
    fn test_get_spectra_by_empty_ids() {
        let storage = SpectrumStorage::new();

        let result = storage.get_spectra_by_ids(&[]);
        assert!(result.is_ok());

        let spectra = result.unwrap();
        assert_eq!(spectra.len(), 0);
    }

    #[test]
    fn test_update_spectrum() {
        let storage = SpectrumStorage::new();
        let spectrum = create_test_spectrum();
        let id = spectrum.id;

        storage.create_spectrum(spectrum).unwrap();

        let baseline = vec![99.5, 100.5, 101.5, 102.5, 103.5, 104.5];
        let corrected = vec![0.5, 0.5, 0.5, 0.5, 0.5, 0.5];

        let result = storage.update_spectrum(id, baseline.clone(), corrected.clone());
        assert!(result.is_ok());

        let updated = result.unwrap();
        assert_eq!(updated.baseline, Some(baseline));
        assert_eq!(updated.corrected, Some(corrected));
    }

    #[test]
    fn test_update_nonexistent_spectrum() {
        let storage = SpectrumStorage::new();
        let id = Uuid::new_v4();

        let baseline = vec![99.5];
        let corrected = vec![0.5];

        let result = storage.update_spectrum(id, baseline, corrected);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_delete_spectrum() {
        let storage = SpectrumStorage::new();
        let spectrum = create_test_spectrum();
        let id = spectrum.id;

        storage.create_spectrum(spectrum).unwrap();

        // Verify it exists
        assert!(storage.get_spectrum(id).is_ok());

        // Delete it
        let result = storage.delete_spectrum(id);
        assert!(result.is_ok());

        // Verify it's gone
        assert!(storage.get_spectrum(id).is_err());
    }

    #[test]
    fn test_delete_nonexistent_spectrum() {
        let storage = SpectrumStorage::new();
        let id = Uuid::new_v4();

        let result = storage.delete_spectrum(id);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_multiple_spectra() {
        let storage = SpectrumStorage::new();

        let spectrum1 = create_test_spectrum();
        let spectrum2 = Spectrum::new(
            "test2.txt".to_string(),
            300,
            305,
            1,
            vec![200, 201, 202, 203, 204, 205],
        );

        let id1 = spectrum1.id;
        let id2 = spectrum2.id;

        storage.create_spectrum(spectrum1).unwrap();
        storage.create_spectrum(spectrum2).unwrap();

        let retrieved1 = storage.get_spectrum(id1).unwrap();
        let retrieved2 = storage.get_spectrum(id2).unwrap();

        assert_eq!(retrieved1.filename, "test.txt");
        assert_eq!(retrieved2.filename, "test2.txt");
        assert_ne!(retrieved1.id, retrieved2.id);
    }
}
