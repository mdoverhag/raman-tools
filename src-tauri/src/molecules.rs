use serde::{Deserialize, Serialize};

// Scientific notation conventions require specific capitalization:
// - Most Raman molecules are fully uppercase acronyms (DTNB, MBA, TFMBA)
// - Target molecules follow standard biochemistry naming:
//   - IgG: Immunoglobulin G (not IGG or Igg)
//   - BSA: Bovine Serum Albumin (fully uppercase)
//   - HER2: Human Epidermal growth factor Receptor 2 (fully uppercase)
//   - EpCAM: Epithelial Cell Adhesion Molecule (mixed case, not EPCAM)
//   - TROP2: Trophoblast cell-surface antigen 2 (fully uppercase)
// We allow the clippy lint to preserve these scientific conventions.

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RamanMolecule {
    DTNB,
    MBA,
    TFMBA,
}

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TargetMolecule {
    IgG,
    BSA,
    HER2,
    EpCAM,
    TROP2,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MoleculePair {
    pub raman: RamanMolecule,
    pub target: TargetMolecule,
}

#[cfg(test)]
impl MoleculePair {
    pub fn new(raman: RamanMolecule, target: TargetMolecule) -> Self {
        Self { raman, target }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_molecule_pair_creation() {
        let pair = MoleculePair::new(RamanMolecule::DTNB, TargetMolecule::HER2);
        assert_eq!(pair.raman, RamanMolecule::DTNB);
        assert_eq!(pair.target, TargetMolecule::HER2);
    }

    #[test]
    fn test_serialization() {
        let pair = MoleculePair::new(RamanMolecule::MBA, TargetMolecule::IgG);
        let json = serde_json::to_string(&pair).unwrap();
        assert!(json.contains("\"raman\":\"MBA\""));
        assert!(json.contains("\"target\":\"IgG\""));

        let deserialized: MoleculePair = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized, pair);

        // Test EpCAM specifically to ensure mixed-case works
        let pair2 = MoleculePair::new(RamanMolecule::DTNB, TargetMolecule::EpCAM);
        let json2 = serde_json::to_string(&pair2).unwrap();
        assert!(json2.contains("\"raman\":\"DTNB\""));
        assert!(json2.contains("\"target\":\"EpCAM\""));
    }
}
