"""
Experiment summary — printed output and JSON export.

Provides a single function that replaces both the old print_experiment_summary
(multiplex only) and the ad-hoc print blocks in singleplex scripts, with
JSON export for version-controlled experiment tracking.
"""

import json
import os
from collections import defaultdict
from datetime import datetime

import numpy as np

from .molecules import get_peak


def experiment_summary(
    experiment,
    samples,
    output_dir,
    references=None,
    deconv_results=None,
):
    """
    Print experiment summary and write a JSON summary file.

    Works for both singleplex (samples only) and multiplex (samples +
    references + deconvolution) experiments.

    Args:
        experiment: Experiment name string. Used for the JSON summary
                    written to summaries/{experiment}.json.
        samples: Dict of sample data (from load_and_process_sample).
        output_dir: Path to the output directory.
        references: Optional dict of reference data (from load_and_process_reference).
        deconv_results: Optional dict of deconvolution results
                        (from normalize_and_deconvolve_samples).
    """
    print("\n" + "=" * 60)
    print("EXPERIMENT COMPLETE")
    print("=" * 60)

    print(f"\nOutput directory: {output_dir}")

    if references:
        print("\nReferences processed:")
        for (molecule, conjugate), data in sorted(references.items()):
            print(f"  {molecule}-{conjugate}: {data['count']} spectra")

    print("\nSamples processed:")
    for sample_key, data in sorted(samples.items(), key=lambda x: x[1]["name"]):
        print(f"  {data['name']}: {data['count']} spectra")

    if deconv_results:
        _print_deconvolution_tables(samples, deconv_results)

    print("\nPlots saved in:")
    print(f"  {output_dir}/")

    summary_path = _summary_path(experiment)
    summary = {
        "results": {
            "references": [],
            "samples": [],
        }
    }

    if references:
        for reference_key, reference_data in sorted(references.items()):
            summary["results"]["references"].append(
                _build_reference_record(reference_key, reference_data)
            )

    for sample_key, sample_data in sorted(samples.items(), key=lambda x: x[1]["name"]):
        sample_deconv = deconv_results.get(sample_key) if deconv_results else None
        summary["results"]["samples"].append(
            _build_sample_record(sample_key, sample_data, sample_deconv)
        )

    _write_json(summary, summary_path)
    print("\nJSON summary written to:")
    print(f"  {summary_path}")


def _print_deconvolution_tables(samples, deconv_results):
    """Print grouped deconvolution contribution tables."""
    # Group samples by their conjugate types
    groups = defaultdict(list)
    for sample_key, sample_data in samples.items():
        conjugate_set = frozenset(
            conj for _, conj in sample_data["molecule_conjugates"]
        )
        groups[conjugate_set].append(sample_key)

    for conjugate_set in sorted(groups.keys(), key=lambda x: sorted(x)):
        sample_keys = groups[conjugate_set]

        first_sample = samples[sample_keys[0]]
        group_mol_conj = sorted(first_sample["molecule_conjugates"])

        conjugate_list = sorted(conjugate_set)
        if len(conjugate_list) == 1:
            group_name = conjugate_list[0]
        else:
            group_name = "/".join(conjugate_list)

        print(f"\n{'=' * 60}")
        print(f"Deconvolution Results - {group_name} Conjugates")
        print(f"{'=' * 60}")

        header = f"{'Sample':<25}"
        for mol, conj in group_mol_conj:
            label = f"{mol}-{conj}"
            header += f" {label:>20}"
        header += f" {'R²':>14}"
        print(f"\n{header}")
        print("-" * (25 + len(group_mol_conj) * 21 + 15))

        for sample_key in sorted(sample_keys, key=lambda k: samples[k]["name"]):
            display_name = samples[sample_key]["name"]
            replicate_results = deconv_results[sample_key]
            sample_mol_conj = samples[sample_key]["molecule_conjugates"]

            contributions_by_mol_conj = {mc: [] for mc in sample_mol_conj}
            r_squared_values = []

            for result in replicate_results:
                for mc in sample_mol_conj:
                    contributions_by_mol_conj[mc].append(
                        result["contributions"][mc]
                    )
                r_squared_values.append(result["metrics"]["r_squared"])

            row = f"{display_name:<25}"
            for mc in group_mol_conj:
                mean = np.mean(contributions_by_mol_conj[mc])
                std = np.std(contributions_by_mol_conj[mc])
                row += f" {f'{mean:.1f}±{std:.1f}%':>20}"

            r2_mean = np.mean(r_squared_values)
            r2_std = np.std(r_squared_values)
            row += f" {f'{r2_mean:.3f}±{r2_std:.3f}':>14}"
            print(row)


def _label_molecule_conjugate(molecule, conjugate):
    """Build a stable label for a molecule-conjugate pair."""
    return f"{molecule}-{conjugate}"


def _summary_path(experiment):
    """Build the summary path for this experiment run."""
    run_id = os.environ.get("RAMAN_SUMMARY_RUN_ID")
    if not run_id:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join("summaries", run_id, f"{experiment}.json")


def _round_float(value, decimals=3):
    """Round a numeric value for stable JSON output."""
    return round(float(value), decimals)


def _distribution_stats(values, decimals=3):
    """Summarize a numeric distribution for diff-friendly JSON."""
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return None

    std = np.std(arr, ddof=1) if arr.size > 1 else 0.0

    return {
        "n": int(arr.size),
        "mean": _round_float(np.mean(arr), decimals),
        "std": _round_float(std, decimals),
        "median": _round_float(np.median(arr), decimals),
        "min": _round_float(np.min(arr), decimals),
        "max": _round_float(np.max(arr), decimals),
        "p05": _round_float(np.percentile(arr, 5), decimals),
        "p95": _round_float(np.percentile(arr, 95), decimals),
    }


def _series_record(values, decimals=3):
    """Build a values + stats record for a numeric series."""
    return {
        "values": [_round_float(value, decimals) for value in values],
        "stats": _distribution_stats(values, decimals),
    }


def _extract_peak_values(replicates, peak_wn):
    """Extract corrected intensities at a peak from replicate spectra."""
    values = []
    for replicate in replicates:
        wavenumbers = np.array(replicate["wavenumbers"], dtype=float)
        corrected = np.array(replicate["corrected"], dtype=float)
        peak_idx = np.argmin(np.abs(wavenumbers - peak_wn))
        values.append(corrected[peak_idx])
    return values


def _build_reference_record(reference_key, reference):
    """Build a JSON-ready summary record for a reference spectrum."""
    molecule, conjugate = reference_key
    peak_wn = get_peak(molecule)
    replicates = reference.get("replicates", [])

    if replicates:
        peak_values = _extract_peak_values(replicates, peak_wn)
    else:
        wavenumbers = np.array(reference["wavenumbers"], dtype=float)
        corrected_avg = np.array(reference["corrected_avg"], dtype=float)
        peak_idx = np.argmin(np.abs(wavenumbers - peak_wn))
        peak_values = [corrected_avg[peak_idx]]

    return {
        "reference_name": _label_molecule_conjugate(molecule, conjugate),
        "molecule": molecule,
        "conjugate": conjugate,
        "count": reference["count"],
        "processing": reference.get("processing", {}),
        "peak_intensity": {
            "wavenumber": peak_wn,
            **_series_record(peak_values),
        },
    }


def _build_sample_record(sample_key, sample, replicate_results=None):
    """Build a JSON-ready summary record for a single sample.

    Args:
        sample_key: Stable sample key from the experiment script.
        sample: Sample dict from load_and_process_sample.
        replicate_results: Optional list of per-replicate deconvolution
                           results for this sample.
    """
    peak_intensities = {}

    for molecule, conjugate in sample["molecule_conjugates"]:
        peak_wn = get_peak(molecule)
        peak_values = _extract_peak_values(sample["replicates"], peak_wn)
        peak_intensities[_label_molecule_conjugate(molecule, conjugate)] = {
            "wavenumber": peak_wn,
            **_series_record(peak_values),
        }

    record = {
        "sample_key": sample_key,
        "sample_name": sample["name"],
        "count": sample["count"],
        "molecule_conjugates": [list(mc) for mc in sample["molecule_conjugates"]],
        "processing": sample.get("processing", {}),
        "peak_intensities": peak_intensities,
    }

    if replicate_results:
        mol_conjs = sample["molecule_conjugates"]
        contributions = {mc: [] for mc in mol_conjs}
        coefficients = {mc: [] for mc in mol_conjs}
        norm_factors = []
        r_squared_values = []
        rmse_values = []
        residual_norm_values = []

        for result in replicate_results:
            for mc in mol_conjs:
                contributions[mc].append(result["contributions"][mc])
                coefficients[mc].append(result["coefficients"][mc])
            if result["norm_factor"] is not None:
                norm_factors.append(result["norm_factor"])
            r_squared_values.append(result["metrics"]["r_squared"])
            rmse_values.append(result["metrics"]["rmse"])
            residual_norm_values.append(result["metrics"]["residual_norm"])

        if norm_factors:
            record["normalization"] = {
                "norm_factors": _series_record(norm_factors),
            }

        record["deconvolution"] = {
            "contributions": {},
            "coefficients": {},
            "fit_metrics": {
                "r_squared": _series_record(r_squared_values),
                "rmse": _series_record(rmse_values),
                "residual_norm": _series_record(residual_norm_values),
            },
        }
        for mol, conj in mol_conjs:
            label = _label_molecule_conjugate(mol, conj)
            record["deconvolution"]["contributions"][label] = _series_record(
                contributions[(mol, conj)]
            )
            record["deconvolution"]["coefficients"][label] = _series_record(
                coefficients[(mol, conj)]
            )

    return record


def _write_json(summary, path):
    """Write summary JSON with indent=2 and trailing newline."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
