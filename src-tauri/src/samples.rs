use crate::molecules::MoleculePair;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Mutex;
use uuid::Uuid;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Sample {
    pub id: Uuid,
    pub name: String,
    pub molecule_pairs: Vec<MoleculePair>,
    pub spectrum_ids: Vec<Uuid>,
    pub average_intensities: Option<Vec<f64>>,
    pub average_corrected: Option<Vec<f64>>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpdateSampleData {
    pub name: Option<String>,
    pub molecule_pairs: Option<Vec<MoleculePair>>,
}

impl Sample {
    pub fn new(name: String) -> Self {
        Self {
            id: Uuid::new_v4(),
            name,
            molecule_pairs: Vec::new(),
            spectrum_ids: Vec::new(),
            average_intensities: None,
            average_corrected: None,
        }
    }
}

// Thread-safe storage for samples
pub struct SampleStorage {
    samples: Mutex<HashMap<Uuid, Sample>>,
}

impl SampleStorage {
    pub fn new() -> Self {
        Self {
            samples: Mutex::new(HashMap::new()),
        }
    }

    pub fn create_sample(&self, name: String) -> Result<Sample, String> {
        let sample = Sample::new(name);
        let mut samples = self.samples.lock().map_err(|e| e.to_string())?;
        samples.insert(sample.id, sample.clone());
        Ok(sample)
    }

    pub fn list_samples(&self) -> Result<Vec<Sample>, String> {
        let samples = self.samples.lock().map_err(|e| e.to_string())?;
        Ok(samples.values().cloned().collect())
    }

    pub fn get_sample(&self, id: Uuid) -> Result<Sample, String> {
        let samples = self.samples.lock().map_err(|e| e.to_string())?;
        samples
            .get(&id)
            .cloned()
            .ok_or_else(|| format!("Sample with id {} not found", id))
    }

    pub fn update_sample(&self, id: Uuid, updates: UpdateSampleData) -> Result<Sample, String> {
        let mut samples = self.samples.lock().map_err(|e| e.to_string())?;

        let sample = samples
            .get_mut(&id)
            .ok_or_else(|| format!("Sample with id {} not found", id))?;

        if let Some(name) = updates.name {
            sample.name = name;
        }
        if let Some(molecule_pairs) = updates.molecule_pairs {
            sample.molecule_pairs = molecule_pairs;
        }

        Ok(sample.clone())
    }

    pub fn delete_sample(&self, id: Uuid) -> Result<(), String> {
        let mut samples = self.samples.lock().map_err(|e| e.to_string())?;
        samples
            .remove(&id)
            .ok_or_else(|| format!("Sample with id {} not found", id))?;
        Ok(())
    }

    pub fn add_spectra_to_sample(&self, id: Uuid, spectrum_ids: Vec<Uuid>) -> Result<(), String> {
        let mut samples = self.samples.lock().map_err(|e| e.to_string())?;
        let sample = samples
            .get_mut(&id)
            .ok_or_else(|| format!("Sample with id {} not found", id))?;

        // Add the new spectrum IDs to the sample
        sample.spectrum_ids.extend(spectrum_ids);

        Ok(())
    }

    pub fn update_average_spectra(
        &self,
        id: Uuid,
        average_intensities: Vec<f64>,
        average_corrected: Vec<f64>,
    ) -> Result<(), String> {
        let mut samples = self.samples.lock().map_err(|e| e.to_string())?;
        let sample = samples
            .get_mut(&id)
            .ok_or_else(|| format!("Sample with id {} not found", id))?;

        sample.average_intensities = Some(average_intensities);
        sample.average_corrected = Some(average_corrected);

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecules::{RamanMolecule, TargetMolecule};

    #[test]
    fn test_create_sample() {
        let storage = SampleStorage::new();
        let result = storage.create_sample("Test Sample".to_string());

        assert!(result.is_ok());
        let sample = result.unwrap();
        assert_eq!(sample.name, "Test Sample");
        assert!(sample.molecule_pairs.is_empty());
        assert!(sample.spectrum_ids.is_empty());
        assert!(sample.average_intensities.is_none());
        assert!(sample.average_corrected.is_none());
    }

    #[test]
    fn test_list_samples() {
        let storage = SampleStorage::new();

        // Initially empty
        let result = storage.list_samples();
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);

        // Add two samples
        storage.create_sample("Sample 1".to_string()).unwrap();
        storage.create_sample("Sample 2".to_string()).unwrap();

        let result = storage.list_samples();
        assert!(result.is_ok());
        let samples = result.unwrap();
        assert_eq!(samples.len(), 2);
    }

    #[test]
    fn test_get_sample() {
        let storage = SampleStorage::new();
        let created = storage.create_sample("Test Sample".to_string()).unwrap();

        // Get existing sample
        let result = storage.get_sample(created.id);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().name, "Test Sample");

        // Get non-existing sample
        let result = storage.get_sample(Uuid::new_v4());
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_update_sample() {
        let storage = SampleStorage::new();
        let created = storage.create_sample("Original Name".to_string()).unwrap();

        // Update name only
        let updates = UpdateSampleData {
            name: Some("Updated Name".to_string()),
            molecule_pairs: None,
        };

        let result = storage.update_sample(created.id, updates);
        assert!(result.is_ok());
        let updated = result.unwrap();
        assert_eq!(updated.name, "Updated Name");

        // Update molecule pairs
        let pairs = vec![
            MoleculePair::new(RamanMolecule::DTNB, TargetMolecule::HER2),
            MoleculePair::new(RamanMolecule::MBA, TargetMolecule::IgG),
        ];
        let updates = UpdateSampleData {
            name: None,
            molecule_pairs: Some(pairs.clone()),
        };

        let result = storage.update_sample(created.id, updates);
        assert!(result.is_ok());
        let updated = result.unwrap();
        assert_eq!(updated.molecule_pairs, pairs);
        assert_eq!(updated.name, "Updated Name"); // Name unchanged
    }

    #[test]
    fn test_update_nonexistent_sample() {
        let storage = SampleStorage::new();
        let updates = UpdateSampleData {
            name: Some("Test".to_string()),
            molecule_pairs: None,
        };

        let result = storage.update_sample(Uuid::new_v4(), updates);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_delete_sample() {
        let storage = SampleStorage::new();
        let created = storage.create_sample("To Delete".to_string()).unwrap();

        // Delete existing sample
        let result = storage.delete_sample(created.id);
        assert!(result.is_ok());

        // Verify it's gone
        let result = storage.get_sample(created.id);
        assert!(result.is_err());

        // Try to delete again
        let result = storage.delete_sample(created.id);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_delete_nonexistent_sample() {
        let storage = SampleStorage::new();
        let result = storage.delete_sample(Uuid::new_v4());
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_add_spectra_to_sample() {
        let storage = SampleStorage::new();
        let sample = storage.create_sample("Test Sample".to_string()).unwrap();

        // Add first batch of spectrum IDs
        let spectrum_ids_1 = vec![Uuid::new_v4(), Uuid::new_v4()];
        let result = storage.add_spectra_to_sample(sample.id, spectrum_ids_1.clone());
        assert!(result.is_ok());

        // Verify the spectra were added
        let updated_sample = storage.get_sample(sample.id).unwrap();
        assert_eq!(updated_sample.spectrum_ids.len(), 2);
        assert_eq!(updated_sample.spectrum_ids, spectrum_ids_1);

        // Add more spectrum IDs (should extend, not replace)
        let spectrum_ids_2 = vec![Uuid::new_v4(), Uuid::new_v4(), Uuid::new_v4()];
        let result = storage.add_spectra_to_sample(sample.id, spectrum_ids_2.clone());
        assert!(result.is_ok());

        // Verify all spectra are present
        let updated_sample = storage.get_sample(sample.id).unwrap();
        assert_eq!(updated_sample.spectrum_ids.len(), 5);

        // Check that original IDs are still there
        assert!(updated_sample.spectrum_ids[0..2]
            .iter()
            .all(|id| spectrum_ids_1.contains(id)));
        assert!(updated_sample.spectrum_ids[2..5]
            .iter()
            .all(|id| spectrum_ids_2.contains(id)));
    }

    #[test]
    fn test_add_spectra_to_nonexistent_sample() {
        let storage = SampleStorage::new();
        let spectrum_ids = vec![Uuid::new_v4()];

        let result = storage.add_spectra_to_sample(Uuid::new_v4(), spectrum_ids);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_add_empty_spectra_list() {
        let storage = SampleStorage::new();
        let sample = storage.create_sample("Test Sample".to_string()).unwrap();

        // Add empty list (should succeed but not change anything)
        let result = storage.add_spectra_to_sample(sample.id, Vec::new());
        assert!(result.is_ok());

        let updated_sample = storage.get_sample(sample.id).unwrap();
        assert_eq!(updated_sample.spectrum_ids.len(), 0);
    }

    #[test]
    fn test_update_average_spectra() {
        let storage = SampleStorage::new();
        let sample = storage.create_sample("Test Sample".to_string()).unwrap();

        // Initially no average spectra
        assert!(sample.average_intensities.is_none());
        assert!(sample.average_corrected.is_none());

        // Update with average spectra
        let avg_intensities = vec![100.5, 101.5, 102.5, 103.5];
        let avg_corrected = vec![95.5, 96.5, 97.5, 98.5];
        let result = storage.update_average_spectra(
            sample.id,
            avg_intensities.clone(),
            avg_corrected.clone(),
        );
        assert!(result.is_ok());

        // Verify they were updated
        let updated_sample = storage.get_sample(sample.id).unwrap();
        assert_eq!(updated_sample.average_intensities, Some(avg_intensities));
        assert_eq!(updated_sample.average_corrected, Some(avg_corrected));
    }

    #[test]
    fn test_update_average_spectra_nonexistent_sample() {
        let storage = SampleStorage::new();
        let avg_intensities = vec![100.5, 101.5];
        let avg_corrected = vec![95.5, 96.5];

        let result = storage.update_average_spectra(Uuid::new_v4(), avg_intensities, avg_corrected);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }
}
