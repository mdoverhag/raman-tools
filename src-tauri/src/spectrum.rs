use serde::{Deserialize, Serialize};
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
    // Not serialized when None to save disk space
    #[serde(skip_serializing_if = "Option::is_none")]
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
