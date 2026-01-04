"""
Raman molecule configuration.

Defines the properties of each Raman-active molecule used in SERS experiments,
including their characteristic peak wavenumbers and plotting colors.
"""

# Configuration for Raman-active molecules
RAMAN_MOLECULES = {
    "MBA": {
        "peak": 1078,  # cm⁻¹
        "color": "blue"
    },
    "DTNB": {
        "peak": 1335,  # cm⁻¹
        "color": "green"
    },
    "TFMBA": {
        "peak": 1377,  # cm⁻¹
        "color": "orange"
    }
}

# Configuration for conjugate colors in histograms
CONJUGATE_COLORS = {
    "BSA": "black",
    "EpCAM": "blue",
    "HER2": "blue",
    "TROP2": "blue",
    "PD-L1": "blue",
    "IgG": "red"
}


def get_peak(molecule: str) -> int:
    """
    Get the characteristic peak wavenumber for a molecule.

    Args:
        molecule: Molecule name (e.g., "MBA", "DTNB", "TFMBA")

    Returns:
        Peak wavenumber in cm⁻¹

    Raises:
        KeyError: If molecule is not defined
    """
    return RAMAN_MOLECULES[molecule]["peak"]


def get_color(molecule: str) -> str:
    """
    Get the plotting color for a molecule.

    Args:
        molecule: Molecule name (e.g., "MBA", "DTNB", "TFMBA")

    Returns:
        Color name for plotting

    Raises:
        KeyError: If molecule is not defined
    """
    return RAMAN_MOLECULES[molecule]["color"]


def get_all_molecules() -> list[str]:
    """
    Get list of all defined Raman molecules.

    Returns:
        List of molecule names
    """
    return list(RAMAN_MOLECULES.keys())
