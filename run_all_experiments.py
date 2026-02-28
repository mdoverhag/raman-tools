#!/usr/bin/env python3
"""Run all experiment scripts in experiments/, using all available cores."""

import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

WORKERS = os.cpu_count()

experiments_dir = Path(__file__).parent / "experiments"
scripts = sorted(experiments_dir.glob("*.py"))

print(f"Found {len(scripts)} experiments, running {WORKERS} at a time\n")


def run_script(script):
    result = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
    )
    return script.name, result


failed = []
with ThreadPoolExecutor(max_workers=WORKERS) as pool:
    futures = {pool.submit(run_script, s): s for s in scripts}
    for i, future in enumerate(as_completed(futures), 1):
        name, result = future.result()
        status = "ok" if result.returncode == 0 else "FAILED"
        print(f"[{i}/{len(scripts)}] {name} — {status}")
        if result.returncode != 0:
            failed.append((name, result))

if failed:
    print(f"\n{len(failed)} experiment(s) failed:\n")
    for name, result in failed:
        print(f"{'=' * 60}")
        print(f"{name} (exit code {result.returncode})")
        print(f"{'=' * 60}")
        if result.stderr:
            print(result.stderr)
    sys.exit(1)

print(f"\nAll {len(scripts)} experiments completed successfully.")
