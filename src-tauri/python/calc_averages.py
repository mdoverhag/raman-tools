#!/usr/bin/env python3
"""
Calculate average and standard deviation of multiple spectra.
"""

import sys
import json
import numpy as np


def calculate_average(spectra):
    """
    Calculate the average and standard deviation of multiple spectra.

    Args:
        spectra: List of spectra, where each spectrum is a list of intensities

    Returns:
        Dictionary with average, standard deviation, and count
    """
    if not spectra:
        raise ValueError("No spectra provided for averaging")

    # Convert to numpy array for efficient computation
    spectra_array = np.array(spectra)

    # Calculate average across spectra (axis 0)
    average = np.mean(spectra_array, axis=0)

    # Calculate standard deviation
    std_dev = np.std(spectra_array, axis=0)

    return {
        "average": average.tolist(),
        "stdDev": std_dev.tolist(),
        "count": len(spectra)
    }


def main():
    # Read input from stdin
    input_data = sys.stdin.read()

    try:
        # Parse JSON input
        data = json.loads(input_data)

        # Extract spectra
        spectra = data.get("spectra", [])

        # Calculate average
        result = calculate_average(spectra)

        # Output JSON result
        print(json.dumps(result))

    except Exception as e:
        # Return error as JSON
        error_result = {
            "error": str(e)
        }
        print(json.dumps(error_result), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()