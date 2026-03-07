#!/usr/bin/env python3
"""Generate HTML reports from summary JSON files or run directories."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Render one summary JSON/run or compare two summary JSON files/run "
            "directories and generate an HTML report."
        )
    )
    parser.add_argument("old_path", help="Summary JSON file or run directory")
    parser.add_argument(
        "new_path",
        nargs="?",
        help=(
            "Optional second summary JSON file or run directory. If omitted, "
            "the script renders a single-summary view."
        ),
    )
    parser.add_argument(
        "--compare-previous",
        action="store_true",
        help=(
            "When only one path is provided, infer the previous comparable "
            "summary and generate a diff instead of a standalone report."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output HTML path. Defaults to "
            "summaries/reports/<name>.html or <old>_vs_<new>.html"
        ),
    )
    return parser.parse_args()


def load_bundle(path: Path) -> dict[str, dict]:
    """Load one or more experiment summaries from a file or directory."""
    if not path.exists():
        raise FileNotFoundError(f"Summary path not found: {path}")

    if path.is_file():
        if path.suffix != ".json":
            raise ValueError(f"Expected a JSON file: {path}")
        with open(path, "r") as f:
            return {path.stem: json.load(f)}

    summaries = {}
    for json_path in sorted(path.glob("*.json")):
        with open(json_path, "r") as f:
            summaries[json_path.stem] = json.load(f)

    if not summaries:
        raise ValueError(f"No JSON summaries found in {path}")

    return summaries


def find_summaries_root(path: Path) -> Path | None:
    """Find the nearest ancestor named 'summaries'."""
    search_path = path if path.is_dir() else path.parent
    for candidate in [search_path, *search_path.parents]:
        if candidate.name == "summaries":
            return candidate
    return None


def default_output_path(old_path: Path, new_path: Path | None) -> Path:
    """Return the default report output path."""
    if new_path is None:
        filename = f"{old_path.stem}.html"
    else:
        filename = f"{old_path.stem}_vs_{new_path.stem}.html"
    summaries_dir = find_summaries_root(old_path) or (
        find_summaries_root(new_path) if new_path is not None else None
    )
    if summaries_dir is None:
        summaries_dir = Path("summaries")
    return summaries_dir / "reports" / filename


def run_directories(parent: Path) -> list[Path]:
    """Return summary run directories under a parent path."""
    dirs = []
    for child in parent.iterdir():
        if not child.is_dir() or child.name in {"diffs", "reports"}:
            continue
        if any(child.glob("*.json")):
            dirs.append(child)
    return sorted(dirs)


def infer_previous_run(run_dir: Path) -> Path:
    """Infer the previous run directory from sibling run directories."""
    parent = run_dir.parent
    candidates = [candidate for candidate in run_directories(parent) if candidate != run_dir]
    older = [candidate for candidate in candidates if candidate.name < run_dir.name]

    if older:
        return older[-1]

    raise ValueError(f"Could not infer previous run directory for {run_dir}")


def infer_previous_file(summary_file: Path) -> Path:
    """Infer the previous summary file with the same name from sibling run dirs."""
    summaries_root = find_summaries_root(summary_file)
    if summaries_root is None:
        raise ValueError(
            f"Cannot infer comparison target for {summary_file}: not under summaries/"
        )

    current_run = summary_file.parent
    older_runs = [
        candidate
        for candidate in run_directories(current_run.parent)
        if candidate != current_run and candidate.name < current_run.name
    ]

    for run_dir in reversed(older_runs):
        candidate = run_dir / summary_file.name
        if candidate.exists():
            return candidate

    raise ValueError(
        f"Could not infer previous summary file for {summary_file.name} from {current_run.parent}"
    )


def resolve_paths(
    old_arg: str, new_arg: str | None, compare_previous: bool
) -> tuple[Path, Path | None]:
    """Resolve input paths, inferring the second path when omitted."""
    first = Path(old_arg)
    if new_arg is not None:
        return first, Path(new_arg)

    if not compare_previous:
        return first, None

    if first.is_dir():
        try:
            return infer_previous_run(first), first
        except ValueError:
            return first, None
    if first.is_file():
        try:
            return infer_previous_file(first), first
        except ValueError:
            return first, None

    raise FileNotFoundError(f"Summary path not found: {first}")


def record_identifier(record_type: str, record: dict) -> str:
    """Return a stable identifier for a sample or reference record."""
    if record_type == "references":
        return record["reference_name"]
    return record.get("sample_key") or record["sample_name"]


def to_record_map(summary: dict, record_type: str) -> dict[str, dict]:
    """Convert summary records to an identifier -> record mapping."""
    records = summary.get("results", {}).get(record_type, [])
    return {record_identifier(record_type, record): record for record in records}


def flatten_value(value, prefix: str = "") -> dict[str, object]:
    """Flatten nested JSON values into a stable path -> value mapping."""
    flat = {}

    if isinstance(value, dict):
        for key, child in sorted(value.items()):
            child_prefix = f"{prefix}.{key}" if prefix else key
            flat.update(flatten_value(child, child_prefix))
        return flat

    if isinstance(value, list):
        flat[prefix] = json.dumps(value, ensure_ascii=True)
        return flat

    flat[prefix] = value
    return flat


def format_value(value) -> str:
    """Format a value for display in the report."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def field_delta(old_value, new_value) -> tuple[str, str]:
    """Return a human-readable delta and a CSS class for the value change."""
    if old_value == new_value:
        return "", "delta-zero"

    if isinstance(old_value, (int, float)) and isinstance(new_value, (int, float)):
        delta = new_value - old_value
        if delta > 0:
            css_class = "delta-up"
        elif delta < 0:
            css_class = "delta-down"
        else:
            css_class = "delta-zero"

        if old_value != 0:
            pct = delta / old_value * 100
            return f"{delta:+.6g} ({pct:+.2f}%)", css_class
        return f"{delta:+.6g}", css_class

    if old_value is None and new_value is not None:
        return "added", "delta-added"
    if old_value is not None and new_value is None:
        return "removed", "delta-removed"
    return "changed", "delta-text"


def build_field_rows(old_record: dict | None, new_record: dict | None) -> list[dict]:
    """Build field-level comparison rows for a record."""
    old_flat = flatten_value(old_record or {})
    new_flat = flatten_value(new_record or {})

    rows = []
    for field in sorted(set(old_flat) | set(new_flat)):
        old_value = old_flat.get(field)
        new_value = new_flat.get(field)

        if old_value is None:
            status = "added"
        elif new_value is None:
            status = "removed"
        elif old_value != new_value:
            status = "changed"
        else:
            status = "unchanged"

        delta, delta_class = field_delta(old_value, new_value)
        rows.append(
            {
                "field": field,
                "old": old_value,
                "new": new_value,
                "status": status,
                "delta": delta,
                "delta_class": delta_class,
            }
        )

    return rows


def count_row_statuses(rows: list[dict]) -> dict[str, int]:
    """Count rows by status for badges and summaries."""
    counts = {"added": 0, "removed": 0, "changed": 0, "unchanged": 0}
    for row in rows:
        counts[row["status"]] += 1
    return counts


def compare_record_group(old_records: dict[str, dict], new_records: dict[str, dict]) -> list[dict]:
    """Compare a record group and return full per-record diffs."""
    diffs = []
    for name in sorted(set(old_records) | set(new_records)):
        old_record = old_records.get(name)
        new_record = new_records.get(name)

        if old_record is None:
            status = "added"
        elif new_record is None:
            status = "removed"
        else:
            status = "changed" if old_record != new_record else "unchanged"

        rows = build_field_rows(old_record, new_record)
        diffs.append(
            {
                "name": name,
                "status": status,
                "rows": rows,
                "counts": count_row_statuses(rows),
            }
        )

    return diffs


def pair_experiments(
    old_bundle: dict[str, dict], new_bundle: dict[str, dict]
) -> tuple[list[dict], list[str], list[str]]:
    """Pair experiments for comparison."""
    old_names = sorted(old_bundle)
    new_names = sorted(new_bundle)

    if len(old_names) == 1 and len(new_names) == 1:
        old_name = old_names[0]
        new_name = new_names[0]
        return (
            [
                {
                    "display_name": (
                        old_name if old_name == new_name else f"{old_name} -> {new_name}"
                    ),
                    "old_name": old_name,
                    "new_name": new_name,
                    "old_summary": old_bundle[old_name],
                    "new_summary": new_bundle[new_name],
                }
            ],
            [],
            [],
        )

    paired = []
    for name in sorted(set(old_names) & set(new_names)):
        paired.append(
            {
                "display_name": name,
                "old_name": name,
                "new_name": name,
                "old_summary": old_bundle[name],
                "new_summary": new_bundle[name],
            }
        )

    added = sorted(set(new_names) - set(old_names))
    removed = sorted(set(old_names) - set(new_names))
    return paired, added, removed


def compare_bundles(old_bundle: dict[str, dict], new_bundle: dict[str, dict]) -> dict:
    """Compare two bundles of experiment summaries."""
    paired, added_experiments, removed_experiments = pair_experiments(
        old_bundle, new_bundle
    )

    experiments = []
    for pair in paired:
        samples = compare_record_group(
            to_record_map(pair["old_summary"], "samples"),
            to_record_map(pair["new_summary"], "samples"),
        )
        references = compare_record_group(
            to_record_map(pair["old_summary"], "references"),
            to_record_map(pair["new_summary"], "references"),
        )
        experiments.append(
            {
                "display_name": pair["display_name"],
                "old_name": pair["old_name"],
                "new_name": pair["new_name"],
                "samples": samples,
                "references": references,
            }
        )

    return {
        "experiments": experiments,
        "added_experiments": added_experiments,
        "removed_experiments": removed_experiments,
    }


def render_name_list(title: str, names: list[str], css_class: str) -> str:
    """Render a titled list of experiment names."""
    if not names:
        return ""

    items = "\n".join(
        f"<li class=\"{css_class}\">{html.escape(name)}</li>" for name in names
    )
    return (
        f"<div class=\"name-list\">"
        f"<h3>{html.escape(title)}</h3>"
        f"<ul>{items}</ul>"
        f"</div>"
    )


def render_record_diff(record: dict, record_type: str) -> str:
    """Render the full field list for one sample or reference."""
    badge_bits = []
    if record["status"] != "unchanged":
        badge_bits.append(record["status"])
    for status in ("changed", "added", "removed"):
        count = record["counts"][status]
        if count:
            badge_bits.append(f"{count} {status}")

    badges = " ".join(
        f"<span class=\"badge badge-{status.split()[1] if ' ' in status else status}\">"
        f"{html.escape(status)}</span>"
        for status in badge_bits
    )

    row_html = []
    for row in record["rows"]:
        row_html.append(
            "<tr>"
            f"<td class=\"field\">{html.escape(row['field'])}</td>"
            f"<td class=\"old old-{row['status']}\"><code>{html.escape(format_value(row['old']))}</code></td>"
            f"<td class=\"new new-{row['status']}\"><code>{html.escape(format_value(row['new']))}</code></td>"
            f"<td class=\"{row['delta_class']}\">{html.escape(row['delta'])}</td>"
            "</tr>"
        )

    status_class = f"record-{record['status']}"
    title = f"{record_type[:-1].title()} {record['name']}"

    return (
        f"<details class=\"record {status_class}\">"
        f"<summary>{html.escape(title)} {badges}</summary>"
        "<table class=\"record-table\">"
        "<thead><tr><th>Field</th><th>Old</th><th>New</th><th>Delta</th></tr></thead>"
        f"<tbody>{''.join(row_html)}</tbody>"
        "</table>"
        "</details>"
    )


def render_single_record(record: dict, record_type: str) -> str:
    """Render one sample or reference without comparison columns."""
    rows = []
    for field, value in flatten_value(record).items():
        rows.append(
            "<tr>"
            f"<td class=\"field\">{html.escape(field)}</td>"
            f"<td class=\"single-value\"><code>{html.escape(format_value(value))}</code></td>"
            "</tr>"
        )

    title = f"{record_type[:-1].title()} {record_identifier(record_type, record)}"
    return (
        "<details class=\"record record-unchanged\">"
        f"<summary>{html.escape(title)}</summary>"
        "<table class=\"record-table single-table\">"
        "<thead><tr><th>Field</th><th>Value</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody>"
        "</table>"
        "</details>"
    )


def render_record_group(records: list[dict], group_name: str) -> str:
    """Render all records in one group."""
    if not records:
        return ""

    changed = sum(1 for record in records if record["status"] == "changed")
    added = sum(1 for record in records if record["status"] == "added")
    removed = sum(1 for record in records if record["status"] == "removed")
    unchanged = sum(1 for record in records if record["status"] == "unchanged")

    header = (
        "<div class=\"group-header\">"
        f"<h3>{html.escape(group_name.title())}</h3>"
        f"<p>{len(records)} records, {changed} changed, {added} added, "
        f"{removed} removed, {unchanged} unchanged</p>"
        "</div>"
    )

    body = "\n".join(render_record_diff(record, group_name) for record in records)
    return f"<section class=\"record-group\">{header}{body}</section>"


def render_single_record_group(records: list[dict], group_name: str) -> str:
    """Render all records in one group without comparison columns."""
    if not records:
        return ""

    header = (
        "<div class=\"group-header\">"
        f"<h3>{html.escape(group_name.title())}</h3>"
        f"<p>{len(records)} records</p>"
        "</div>"
    )

    body = "\n".join(
        render_single_record(record, group_name) for record in records
    )
    return f"<section class=\"record-group\">{header}{body}</section>"


def render_experiment(experiment: dict) -> str:
    """Render one full experiment diff section."""
    title = experiment["display_name"]
    if experiment["old_name"] != experiment["new_name"]:
        subtitle = (
            f"<p class=\"experiment-pair\">Old name: {html.escape(experiment['old_name'])} "
            f"| New name: {html.escape(experiment['new_name'])}</p>"
        )
    else:
        subtitle = ""

    return (
        "<details class=\"experiment\">"
        f"<summary><span>{html.escape(title)}</span></summary>"
        "<div class=\"experiment-body\">"
        f"{subtitle}"
        f"{render_record_group(experiment['references'], 'references')}"
        f"{render_record_group(experiment['samples'], 'samples')}"
        "</div>"
        "</details>"
    )


def render_single_experiment(name: str, summary: dict) -> str:
    """Render one full experiment section without comparison."""
    return (
        "<details class=\"experiment\">"
        f"<summary><span>{html.escape(name)}</span></summary>"
        "<div class=\"experiment-body\">"
        f"{render_single_record_group(summary.get('results', {}).get('references', []), 'references')}"
        f"{render_single_record_group(summary.get('results', {}).get('samples', []), 'samples')}"
        "</div>"
        "</details>"
    )


def diff_styles() -> str:
    """Return shared CSS for diff reports."""
    return """
    :root {
      --bg: #f4f1ea;
      --panel: #fffdf8;
      --ink: #1f1a17;
      --muted: #70675f;
      --line: #ddd2c5;
      --up: #1f7a4d;
      --down: #b44332;
      --changed: #a66a00;
      --added-bg: #e5f5e8;
      --removed-bg: #fae0db;
      --changed-bg: #fff1d6;
      --unchanged-bg: #f8f5ef;
      --badge-bg: #eee5d8;
      --code-bg: #f7f2ea;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 28px;
      color: var(--ink);
      font-family: Georgia, "Times New Roman", serif;
      background:
        radial-gradient(circle at top left, #ebe6dc 0%, transparent 28%),
        linear-gradient(180deg, #f7f4ef 0%, var(--bg) 100%);
    }
    main {
      max-width: 1400px;
      margin: 0 auto;
    }
    h1, h2, h3 {
      margin-top: 0;
    }
    .intro, .experiment, .record-group {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 14px;
      box-shadow: 0 10px 24px rgba(31, 26, 23, 0.05);
    }
    .intro {
      padding: 20px 22px;
      margin-bottom: 22px;
    }
    .experiment {
      margin-bottom: 22px;
      overflow: hidden;
    }
    .experiment > summary {
      cursor: pointer;
      list-style: none;
      padding: 16px 20px;
      font-size: 1.25rem;
      font-weight: bold;
      background: #f3ece1;
    }
    .experiment > summary::-webkit-details-marker {
      display: none;
    }
    .experiment-body {
      padding: 20px 22px;
    }
    .record-group {
      padding: 16px 18px;
      margin: 16px 0;
    }
    .paths {
      color: var(--muted);
      line-height: 1.5;
    }
    .name-list ul {
      margin: 8px 0 0;
      padding-left: 20px;
    }
    .added {
      color: var(--up);
    }
    .removed {
      color: var(--down);
    }
    .group-header p, .experiment-pair {
      color: var(--muted);
      margin-top: -6px;
    }
    details.record {
      margin: 14px 0;
      border: 1px solid var(--line);
      border-radius: 12px;
      overflow: hidden;
    }
    details.record > summary {
      cursor: pointer;
      padding: 12px 14px;
      list-style: none;
      font-weight: bold;
    }
    details.record > summary::-webkit-details-marker {
      display: none;
    }
    .record-added > summary {
      background: var(--added-bg);
    }
    .record-removed > summary {
      background: var(--removed-bg);
    }
    .record-changed > summary {
      background: var(--changed-bg);
    }
    .record-unchanged > summary {
      background: var(--unchanged-bg);
    }
    .badge {
      display: inline-block;
      margin-left: 8px;
      padding: 3px 9px;
      border-radius: 999px;
      background: var(--badge-bg);
      color: var(--muted);
      font-size: 0.82rem;
      font-weight: normal;
    }
    .badge-added {
      color: var(--up);
    }
    .badge-removed {
      color: var(--down);
    }
    .badge-changed {
      color: var(--changed);
    }
    table {
      width: 100%;
      border-collapse: collapse;
      table-layout: fixed;
    }
    th, td {
      border-top: 1px solid var(--line);
      padding: 9px 10px;
      text-align: left;
      vertical-align: top;
    }
    th {
      background: #f2eadf;
    }
    td.field {
      width: 28%;
      font-weight: bold;
    }
    td.old, td.new {
      width: 29%;
    }
    td code {
      display: block;
      white-space: pre-wrap;
      word-break: break-word;
      background: var(--code-bg);
      padding: 6px 8px;
      border-radius: 8px;
      font-family: "SFMono-Regular", Menlo, monospace;
      font-size: 0.84rem;
    }
    .old-added code, .new-added code {
      background: var(--added-bg);
    }
    .old-removed code, .new-removed code {
      background: var(--removed-bg);
    }
    .old-changed code, .new-changed code {
      background: var(--changed-bg);
    }
    .delta-up {
      color: var(--up);
      font-weight: bold;
    }
    .delta-down {
      color: var(--down);
      font-weight: bold;
    }
    .delta-added {
      color: var(--up);
      font-weight: bold;
    }
    .delta-removed {
      color: var(--down);
      font-weight: bold;
    }
    .delta-zero, .delta-text {
      color: var(--muted);
    }
    """


def single_styles() -> str:
    """Return shared CSS for single-summary reports."""
    return """
    :root {
      --bg: #f4f1ea;
      --panel: #fffdf8;
      --ink: #1f1a17;
      --muted: #70675f;
      --line: #ddd2c5;
      --code-bg: #f7f2ea;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 28px;
      color: var(--ink);
      font-family: Georgia, "Times New Roman", serif;
      background:
        radial-gradient(circle at top left, #ebe6dc 0%, transparent 28%),
        linear-gradient(180deg, #f7f4ef 0%, var(--bg) 100%);
    }
    main {
      max-width: 1400px;
      margin: 0 auto;
    }
    h1, h2, h3 {
      margin-top: 0;
    }
    .intro, .experiment, .record-group {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 14px;
      box-shadow: 0 10px 24px rgba(31, 26, 23, 0.05);
    }
    .intro {
      padding: 20px 22px;
      margin-bottom: 22px;
    }
    .experiment {
      margin-bottom: 22px;
      overflow: hidden;
    }
    .experiment > summary {
      cursor: pointer;
      list-style: none;
      padding: 16px 20px;
      font-size: 1.25rem;
      font-weight: bold;
      background: #f3ece1;
    }
    .experiment > summary::-webkit-details-marker {
      display: none;
    }
    .experiment-body {
      padding: 20px 22px;
    }
    .record-group {
      padding: 16px 18px;
      margin: 16px 0;
    }
    .paths, .group-header p {
      color: var(--muted);
      line-height: 1.5;
    }
    details.record {
      margin: 14px 0;
      border: 1px solid var(--line);
      border-radius: 12px;
      overflow: hidden;
    }
    details.record > summary {
      cursor: pointer;
      padding: 12px 14px;
      list-style: none;
      font-weight: bold;
      background: #f8f5ef;
    }
    details.record > summary::-webkit-details-marker {
      display: none;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      table-layout: fixed;
    }
    th, td {
      border-top: 1px solid var(--line);
      padding: 9px 10px;
      text-align: left;
      vertical-align: top;
    }
    th {
      background: #f2eadf;
    }
    td.field {
      width: 35%;
      font-weight: bold;
    }
    td.single-value {
      width: 65%;
    }
    td code {
      display: block;
      white-space: pre-wrap;
      word-break: break-word;
      background: var(--code-bg);
      padding: 6px 8px;
      border-radius: 8px;
      font-family: "SFMono-Regular", Menlo, monospace;
      font-size: 0.84rem;
    }
    """


def render_html(old_path: Path, new_path: Path, comparison: dict) -> str:
    """Render the full HTML diff report."""
    experiment_sections = "\n".join(
        render_experiment(experiment) for experiment in comparison["experiments"]
    )
    added_experiments = render_name_list(
        "Added Experiments", comparison["added_experiments"], "added"
    )
    removed_experiments = render_name_list(
        "Removed Experiments", comparison["removed_experiments"], "removed"
    )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Summary Diff Report</title>
  <style>{diff_styles()}</style>
</head>
<body>
  <main>
    <section class="intro">
      <h1>Summary Diff Report</h1>
      <p class="paths">
        <strong>Old:</strong> {html.escape(str(old_path))}<br>
        <strong>New:</strong> {html.escape(str(new_path))}
      </p>
      {added_experiments}
      {removed_experiments}
    </section>
    {experiment_sections}
  </main>
</body>
</html>
"""


def render_single_html(path: Path, bundle: dict[str, dict]) -> str:
    """Render the full HTML report for a single bundle."""
    experiment_sections = "\n".join(
        render_single_experiment(name, summary)
        for name, summary in sorted(bundle.items())
    )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Summary Report</title>
  <style>{single_styles()}</style>
</head>
<body>
  <main>
    <section class="intro">
      <h1>Summary Report</h1>
      <p class="paths"><strong>Source:</strong> {html.escape(str(path))}</p>
    </section>
    {experiment_sections}
  </main>
</body>
</html>
"""


def main() -> int:
    """Entry point."""
    args = parse_args()
    old_path, new_path = resolve_paths(
        args.old_path, args.new_path, args.compare_previous
    )

    if new_path is None:
        bundle = load_bundle(old_path)
        output_path = Path(args.output) if args.output else default_output_path(
            old_path, None
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(render_single_html(old_path, bundle))
        print(f"Report written to: {output_path}")
        return 0

    old_bundle = load_bundle(old_path)
    new_bundle = load_bundle(new_path)
    comparison = compare_bundles(old_bundle, new_bundle)

    output_path = Path(args.output) if args.output else default_output_path(
        old_path, new_path
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_html(old_path, new_path, comparison))

    print(f"Diff report written to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
