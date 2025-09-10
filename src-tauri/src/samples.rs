use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Mutex;
use uuid::Uuid;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Sample {
    pub id: Uuid,
    pub name: String,
    pub raman_molecules: Vec<String>,
    pub target_molecules: Vec<String>,
    pub spectrum_ids: Vec<Uuid>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpdateSampleData {
    pub name: Option<String>,
    pub raman_molecules: Option<Vec<String>>,
    pub target_molecules: Option<Vec<String>>,
}

impl Sample {
    pub fn new(name: String) -> Self {
        Self {
            id: Uuid::new_v4(),
            name,
            raman_molecules: Vec::new(),
            target_molecules: Vec::new(),
            spectrum_ids: Vec::new(),
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
        if let Some(raman_molecules) = updates.raman_molecules {
            sample.raman_molecules = raman_molecules;
        }
        if let Some(target_molecules) = updates.target_molecules {
            sample.target_molecules = target_molecules;
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_sample() {
        let storage = SampleStorage::new();
        let result = storage.create_sample("Test Sample".to_string());

        assert!(result.is_ok());
        let sample = result.unwrap();
        assert_eq!(sample.name, "Test Sample");
        assert!(sample.raman_molecules.is_empty());
        assert!(sample.target_molecules.is_empty());
        assert!(sample.spectrum_ids.is_empty());
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
            raman_molecules: None,
            target_molecules: None,
        };

        let result = storage.update_sample(created.id, updates);
        assert!(result.is_ok());
        let updated = result.unwrap();
        assert_eq!(updated.name, "Updated Name");

        // Update molecules
        let updates = UpdateSampleData {
            name: None,
            raman_molecules: Some(vec!["DTNB".to_string(), "MBA".to_string()]),
            target_molecules: Some(vec!["IgG".to_string()]),
        };

        let result = storage.update_sample(created.id, updates);
        assert!(result.is_ok());
        let updated = result.unwrap();
        assert_eq!(updated.raman_molecules, vec!["DTNB", "MBA"]);
        assert_eq!(updated.target_molecules, vec!["IgG"]);
        assert_eq!(updated.name, "Updated Name"); // Name unchanged
    }

    #[test]
    fn test_update_nonexistent_sample() {
        let storage = SampleStorage::new();
        let updates = UpdateSampleData {
            name: Some("Test".to_string()),
            raman_molecules: None,
            target_molecules: None,
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
}
