#!/usr/bin/env python3
"""
Batch processor for baseline correction with streaming output.
Processes multiple spectra and streams results back as they complete.
"""

import json
import sys
from baseline_correction import process_spectrum


def process_batch():
    """
    Process a batch of spectra, streaming results as JSON lines.

    Input format (via stdin):
    {
        "spectra": [
            {"index": 0, "intensities": [...]},
            {"index": 1, "intensities": [...]},
            ...
        ],
        "params": {
            "denoise": bool,
            "window_size": int,
            "lambda_param": float,
            "p": float,
            "d": int
        }
    }

    Output format (via stdout, one JSON object per line):
    {"type": "progress", "current": 1, "total": 150}
    {"type": "result", "index": 0, "baseline": [...], "corrected": [...], "denoised": [...]}
    {"type": "progress", "current": 2, "total": 150}
    {"type": "result", "index": 1, "baseline": [...], "corrected": [...], "denoised": [...]}
    ...
    {"type": "complete", "processed": 150}
    """
    try:
        # Read the entire batch request
        input_data = json.loads(sys.stdin.read())
        spectra = input_data["spectra"]
        params = input_data.get("params", {})

        total = len(spectra)
        processed = 0

        # Process each spectrum and stream results
        for spectrum_data in spectra:
            index = spectrum_data["index"]
            intensities = spectrum_data["intensities"]

            try:
                # Send progress update
                progress = {
                    "type": "progress",
                    "current": processed + 1,
                    "total": total,
                }
                print(json.dumps(progress), flush=True)

                # Process the spectrum
                corrected, baseline, denoised = process_spectrum(
                    spectrum=intensities,
                    denoise=params.get("denoise", False),
                    window_size=params.get("window_size", 5),
                    lambda_param=params.get("lambda_param", 1e7),
                    p=params.get("p", 0.01),
                    d=params.get("d", 2),
                )

                # Send result
                result = {
                    "type": "result",
                    "index": index,
                    "baseline": baseline,
                    "corrected": corrected,
                }
                if denoised is not None:
                    result["denoised"] = denoised

                print(json.dumps(result), flush=True)
                processed += 1

            except Exception as e:
                # Send error for this specific spectrum but continue with others
                error = {
                    "type": "error",
                    "index": index,
                    "error": str(e),
                    "error_type": type(e).__name__,
                }
                print(json.dumps(error), flush=True)

        # Send completion message
        complete = {"type": "complete", "processed": processed}
        print(json.dumps(complete), flush=True)

    except Exception as e:
        # Fatal error - couldn't parse input or other critical failure
        error = {"type": "fatal_error", "error": str(e), "error_type": type(e).__name__}
        print(json.dumps(error), flush=True)
        sys.exit(1)


if __name__ == "__main__":
    process_batch()
