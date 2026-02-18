"""
Input/Output functions for Raman spectroscopy data.

Handles loading spectrum files from disk and saving plots.
"""

from pathlib import Path


def load_spectrum(filepath: str) -> dict:
    """
    Load a single spectrum file.

    Expected format: Tab-delimited text file with two columns:
    - Column 1: Wavenumber (cm⁻¹)
    - Column 2: Intensity

    Args:
        filepath: Path to the spectrum .txt file

    Returns:
        Dictionary with 'wavenumbers' and 'intensities' as lists

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    path = Path(filepath)

    if not path.exists():
        raise FileNotFoundError(f"Spectrum file not found: {path}")

    wavenumbers = []
    intensities = []

    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Parse tab-delimited data
            parts = line.split("\t")

            if len(parts) != 2:
                raise ValueError(
                    f"Invalid format in {path.name} at line {line_num}: "
                    f"expected 2 columns, got {len(parts)}"
                )

            try:
                wavenumber = float(parts[0])
                intensity = float(parts[1])
            except ValueError as e:
                raise ValueError(f"Invalid data in {path.name} at line {line_num}: {e}")

            wavenumbers.append(wavenumber)
            intensities.append(intensity)

    if not wavenumbers:
        raise ValueError(f"No data found in {path}")

    return {
        "wavenumbers": wavenumbers,
        "intensities": intensities,
        "filename": path.name,
    }


def load_spectra(directory: str) -> list[dict]:
    """
    Load all spectrum files from a directory.

    Scans for all .txt files in the directory and loads them as spectra.

    Args:
        directory: Path to directory containing spectrum .txt files

    Returns:
        List of spectrum dictionaries, sorted by filename

    Raises:
        FileNotFoundError: If directory doesn't exist
        ValueError: If no .txt files found or files have invalid format
    """
    dir_path = Path(directory)

    if not dir_path.exists():
        raise FileNotFoundError(f"Directory not found: {dir_path}")

    if not dir_path.is_dir():
        raise ValueError(f"Not a directory: {dir_path}")

    # Find all .txt files
    txt_files = sorted(dir_path.glob("*.txt"))

    if not txt_files:
        raise ValueError(f"No .txt files found in {directory}")

    spectra = []
    errors = []

    for filepath in txt_files:
        try:
            spectrum = load_spectrum(str(filepath))
            spectra.append(spectrum)
        except Exception as e:
            errors.append(f"{filepath.name}: {e}")

    if errors:
        error_msg = "Errors loading spectra:\n" + "\n".join(errors)
        raise ValueError(error_msg)

    if not spectra:
        raise ValueError(f"No valid spectra loaded from {directory}")

    return spectra


def load_multicolumn_spectra(filepath: str) -> list[dict]:
    """
    Load spectra from a multi-column file.

    Expected format: Tab-delimited text file where:
    - Column 1: Wavenumber (cm⁻¹)
    - Columns 2+: Intensity values (one column per replicate)

    Args:
        filepath: Path to the multi-column spectrum .txt file

    Returns:
        List of spectrum dictionaries (one per column), sorted by column index.
        Each dict has 'wavenumbers', 'intensities', and 'filename' keys.

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid or no data found
    """
    path = Path(filepath)

    if not path.exists():
        raise FileNotFoundError(f"Spectrum file not found: {path}")

    wavenumbers = []
    columns = []  # list of lists, one per replicate

    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(
                    f"Invalid format in {path.name} at line {line_num}: "
                    f"expected at least 3 columns (wavenumber + 2 intensities), "
                    f"got {len(parts)}"
                )

            try:
                wavenumber = float(parts[0])
            except ValueError as e:
                raise ValueError(
                    f"Invalid wavenumber in {path.name} at line {line_num}: {e}"
                )

            # Initialize columns on first data line
            if not columns:
                columns = [[] for _ in range(len(parts) - 1)]
            elif len(parts) - 1 != len(columns):
                raise ValueError(
                    f"Inconsistent column count in {path.name} at line {line_num}: "
                    f"expected {len(columns) + 1} columns, got {len(parts)}"
                )

            wavenumbers.append(wavenumber)

            for i, val in enumerate(parts[1:]):
                try:
                    columns[i].append(float(val))
                except ValueError as e:
                    raise ValueError(
                        f"Invalid intensity in {path.name} at line {line_num}, "
                        f"column {i + 2}: {e}"
                    )

    if not wavenumbers:
        raise ValueError(f"No data found in {path}")

    # Build one spectrum dict per column
    spectra = []
    for i, intensities in enumerate(columns):
        spectra.append({
            "wavenumbers": list(wavenumbers),
            "intensities": intensities,
            "filename": f"{path.stem}_{i + 1:03d}",
        })

    return spectra


def create_output_dir(name: str, base_dir: str = ".") -> str:
    """
    Create an output directory with automatic versioning.

    If the directory already exists, appends a number (name-2, name-3, etc.)

    Args:
        name: Base name for the directory
        base_dir: Parent directory (default: current directory)

    Returns:
        Path to the created directory

    Example:
        >>> create_output_dir("experiment-1")
        "experiment-1"  # created
        >>> create_output_dir("experiment-1")
        "experiment-1-2"  # since experiment-1 exists
    """
    base_path = Path(base_dir)
    base_path.mkdir(parents=True, exist_ok=True)

    # Try the base name first
    output_path = base_path / name

    if not output_path.exists():
        output_path.mkdir(parents=True)
        return str(output_path)

    # If exists, try with incrementing numbers
    counter = 2
    while True:
        versioned_name = f"{name}-{counter}"
        output_path = base_path / versioned_name

        if not output_path.exists():
            output_path.mkdir(parents=True)
            return str(output_path)

        counter += 1

        # Safety check to avoid infinite loops
        if counter > 1000:
            raise RuntimeError(
                f"Could not create output directory: too many versions of '{name}'"
            )


def ensure_output_subdir(output_dir: str, subdir_name: str) -> str:
    """
    Ensure a subdirectory exists within the output directory.

    Args:
        output_dir: Base output directory path
        subdir_name: Name of subdirectory to create

    Returns:
        Path to the subdirectory
    """
    subdir_path = Path(output_dir) / subdir_name
    subdir_path.mkdir(parents=True, exist_ok=True)
    return str(subdir_path)
