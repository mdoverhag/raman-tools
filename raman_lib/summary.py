"""
Experiment summary — printed output and optional JSON export.

Provides a single function that replaces both the old print_experiment_summary
(multiplex only) and the ad-hoc print blocks in singleplex scripts, with
optional JSON export for version-controlled experiment tracking.
"""

import json
import os
from collections import defaultdict

import numpy as np

from .molecules import get_peak


def experiment_summary(
    samples,
    output_dir,
    references=None,
    deconv_results=None,
    summary_path=None,
):
    """
    Print experiment summary and optionally write a JSON summary file.

    Works for both singleplex (samples only) and multiplex (samples +
    references + deconvolution) experiments.

    Args:
        samples: Dict of sample data (from load_and_process_sample).
        output_dir: Path to the output directory.
        references: Optional dict of reference data (from load_and_process_reference).
        deconv_results: Optional dict of deconvolution results
                        (from normalize_and_deconvolve_samples).
        summary_path: If provided, write a JSON summary to this path.
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

    if summary_path:
        records = []
        for sample_key, sample_data in samples.items():
            sample_deconv = deconv_results.get(sample_key) if deconv_results else None
            records.append(_build_sample_record(sample_data, sample_deconv))
        _write_json(records, summary_path)
        print(f"\nJSON summary written to:")
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


def _build_sample_record(sample, replicate_results=None):
    """Build a JSON-ready summary record for a single sample.

    Args:
        sample: Sample dict from load_and_process_sample.
        replicate_results: Optional list of per-replicate deconvolution
                           results for this sample.
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

    record = {
        "sample_name": sample["name"],
        "count": sample["count"],
        "molecule_conjugates": [list(mc) for mc in sample["molecule_conjugates"]],
        "peak_intensities": peak_intensities,
    }

    if replicate_results:
        mol_conjs = sample["molecule_conjugates"]
        contributions = {mc: [] for mc in mol_conjs}
        r_squared_values = []

        for result in replicate_results:
            for mc in mol_conjs:
                contributions[mc].append(result["contributions"][mc])
            r_squared_values.append(result["metrics"]["r_squared"])

        record["deconvolution"] = {}
        for mol, conj in mol_conjs:
            vals = np.array(contributions[(mol, conj)])
            record["deconvolution"][f"{mol}-{conj}"] = {
                "mean": round(float(np.mean(vals)), 1),
                "std": round(float(np.std(vals, ddof=1)), 1),
            }
        r2 = np.array(r_squared_values)
        record["deconvolution"]["r_squared"] = {
            "mean": round(float(np.mean(r2)), 3),
            "std": round(float(np.std(r2, ddof=1)), 3),
        }

    return record


def _write_json(records, path):
    """Write records as JSON with indent=2 and trailing newline."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(records, f, indent=2)
        f.write("\n")
