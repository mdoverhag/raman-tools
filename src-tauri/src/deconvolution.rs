use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Mutex;
use uuid::Uuid;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DeconvolutionRun {
    pub id: Uuid,
    pub multiplex_sample_id: Uuid,
    pub reference_sample_ids: Vec<Uuid>, // Singleplex sample IDs for each Raman molecule
    pub wavenumber_range: (f64, f64),    // e.g., (1000.0, 1500.0)
    pub normalized_multiplex: Option<Vec<f64>>,
    pub normalized_references: Option<HashMap<String, Vec<f64>>>, // Keyed by Raman molecule name
    pub deconvolution_results: Option<DeconvolutionResults>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DeconvolutionResults {
    pub contributions: HashMap<String, f64>, // Raman molecule -> percentage contribution
    pub reconstructed_spectrum: Vec<f64>,
    pub residual: Vec<f64>,
    pub rmse: f64,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct CreateDeconvolutionData {
    pub multiplex_sample_id: Uuid,
    pub reference_sample_ids: Vec<Uuid>,
    pub wavenumber_range: Option<(f64, f64)>, // Default to (1000, 1500) if not provided
}

impl DeconvolutionRun {
    pub fn new(
        multiplex_sample_id: Uuid,
        reference_sample_ids: Vec<Uuid>,
        wavenumber_range: (f64, f64),
    ) -> Self {
        Self {
            id: Uuid::new_v4(),
            multiplex_sample_id,
            reference_sample_ids,
            wavenumber_range,
            normalized_multiplex: None,
            normalized_references: None,
            deconvolution_results: None,
        }
    }
}

pub struct DeconvolutionStorage {
    runs: Mutex<HashMap<Uuid, DeconvolutionRun>>,
}

impl DeconvolutionStorage {
    pub fn new() -> Self {
        Self {
            runs: Mutex::new(HashMap::new()),
        }
    }

    pub fn create_run(
        &self,
        multiplex_sample_id: Uuid,
        reference_sample_ids: Vec<Uuid>,
        wavenumber_range: Option<(f64, f64)>,
    ) -> Result<DeconvolutionRun, String> {
        let range = wavenumber_range.unwrap_or((1000.0, 1500.0));
        let run = DeconvolutionRun::new(multiplex_sample_id, reference_sample_ids, range);

        let mut runs = self.runs.lock().map_err(|e| e.to_string())?;
        runs.insert(run.id, run.clone());

        Ok(run)
    }

    pub fn get_run(&self, id: &Uuid) -> Result<DeconvolutionRun, String> {
        let runs = self.runs.lock().map_err(|e| e.to_string())?;
        runs.get(id)
            .cloned()
            .ok_or_else(|| "Deconvolution run not found".to_string())
    }

    pub fn update_normalized_data(
        &self,
        id: &Uuid,
        normalized_multiplex: Vec<f64>,
        normalized_references: HashMap<String, Vec<f64>>,
    ) -> Result<(), String> {
        let mut runs = self.runs.lock().map_err(|e| e.to_string())?;
        let run = runs
            .get_mut(id)
            .ok_or_else(|| "Deconvolution run not found".to_string())?;

        run.normalized_multiplex = Some(normalized_multiplex);
        run.normalized_references = Some(normalized_references);

        Ok(())
    }
}
