#!/usr/bin/env python3
"""
Calculate average and standard deviation of multiple spectra.
"""

import sys
import json
import numpy as np


def calculate_averages(raw_spectra, corrected_spectra):
    """
    Calculate the average and standard deviation of raw and corrected spectra.

    Args:
        raw_spectra: List of raw intensity spectra
        corrected_spectra: List of baseline-corrected spectra

    Returns:
        Dictionary with averages, standard deviations, and count for both raw and corrected
    """
    if not raw_spectra:
        raise ValueError("No spectra provided for averaging")

    result = {}

    # Calculate raw intensities average
    raw_array = np.array(raw_spectra)
    result["averageIntensities"] = np.mean(raw_array, axis=0).tolist()
    result["stdDevIntensities"] = np.std(raw_array, axis=0).tolist()

    # Calculate corrected average if provided
    if corrected_spectra:
        corrected_array = np.array(corrected_spectra)
        result["averageCorrected"] = np.mean(corrected_array, axis=0).tolist()
        result["stdDevCorrected"] = np.std(corrected_array, axis=0).tolist()
    else:
        result["averageCorrected"] = None
        result["stdDevCorrected"] = None

    result["count"] = len(raw_spectra)

    return result


def main():
    # Read input from stdin
    input_data = sys.stdin.read()

    try:
        # Parse JSON input
        data = json.loads(input_data)

        # Extract raw and corrected spectra
        raw_spectra = data.get("rawSpectra", [])
        corrected_spectra = data.get("correctedSpectra", [])

        # Calculate averages
        result = calculate_averages(raw_spectra, corrected_spectra)

        # Output JSON result
        print(json.dumps(result))

    except Exception as e:
        # Return error as JSON
        error_result = {"error": str(e)}
        print(json.dumps(error_result), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
