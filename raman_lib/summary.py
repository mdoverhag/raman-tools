"""
JSON summary export for experiment tracking.

Produces structured records from sample data so that algorithm changes
(baseline correction, normalization, deconvolution) surface clearly in
version-controlled JSON diffs.
"""

import json
import os

import numpy as np

from .molecules import get_peak


def build_sample_record(sample, experiment):
    """
    Build a JSON-ready summary record for a single sample.

    For each molecule in the sample's molecule_conjugates, extracts the
    corrected intensity at the characteristic peak wavenumber across all
    replicates and computes mean and std.

    Args:
        sample: Sample dict from load_and_process_sample, must contain
                'name', 'count', 'molecule_conjugates', and 'replicates'.
        experiment: Experiment name string (typically the script filename
                    without extension).

    Returns:
        Dict with keys: action, experiment, sample_name, count,
        molecule_conjugates, peak_intensities.
    """
    peak_intensities = {}

    for molecule, conjugate in sample["molecule_conjugates"]:
        peak_wn = get_peak(molecule)
        intensities = []

        for replicate in sample["replicates"]:
            wavenumbers = np.array(replicate["wavenumbers"])
            corrected = np.array(replicate["corrected"])
            peak_idx = np.argmin(np.abs(wavenumbers - peak_wn))
            intensities.append(corrected[peak_idx])

        arr = np.array(intensities)
        peak_intensities[molecule] = {
            "mean": round(float(np.mean(arr)), 1),
            "std": round(float(np.std(arr, ddof=1)), 1),
            "wavenumber": peak_wn,
        }

    return {
        "action": "load_sample",
        "experiment": experiment,
        "sample_name": sample["name"],
        "count": sample["count"],
        "molecule_conjugates": [
            list(mc) for mc in sample["molecule_conjugates"]
        ],
        "peak_intensities": peak_intensities,
    }


def write_summary(records, path):
    """
    Write a list of summary records to a JSON file.

    Creates parent directories if needed. Writes with indent=2 and a
    trailing newline for clean git diffs.

    Args:
        records: List of record dicts (from build_sample_record).
        path: File path to write to.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(records, f, indent=2)
        f.write("\n")
