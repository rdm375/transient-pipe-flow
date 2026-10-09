#!/usr/bin/env python3
"""M10.1 reproducible end-to-end CPU performance sweep."""
import argparse
import csv
import json
import os
import platform
import statistics
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXECUTABLE = ROOT / "build" / "m10_cpu_baseline"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path,
                        default=EXECUTABLE)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--duration", type=float, default=1200.0)
    parser.add_argument("--meshes", nargs="+", type=int,
                        default=[10, 20, 40, 80, 100])
    parser.add_argument("--timesteps", nargs="+", type=float,
                        default=[60.0, 30.0, 15.0])
    parser.add_argument("--output", type=Path,
                        default=ROOT / "benchmarks/results-m10")
    args = parser.parse_args()

    if args.repetitions < 1:
        parser.error("repetitions must be positive")
    executable = args.executable.resolve()
    if not executable.is_file():
        parser.error(f"Executable is missing: {executable}")

    args.output.mkdir(parents=True, exist_ok=True)
    rows = []

    for n in args.meshes:
        for dt in args.timesteps:
            measurements = []
            for repetition in range(args.repetitions):
                command = [str(executable), str(n), str(dt),
                           str(args.duration)]
                started = time.perf_counter()
                result = subprocess.run(
                    command, cwd=ROOT, text=True,
                    capture_output=True, check=True)
                wall = time.perf_counter() - started

                lines = result.stdout.strip().splitlines()
                if len(lines) != 2:
                    raise RuntimeError(
                        f"Unexpected benchmark output: {result.stdout}")
                parsed = next(csv.DictReader(lines))
                row = {
                    key: (int(value) if key in
                          ("n", "steps", "newton_total")
                          else float(value))
                    for key, value in parsed.items()
                }
                row["repetition"] = repetition
                row["wall_seconds"] = wall
                rows.append(row)
                measurements.append(row["cpu_seconds"])

            print(f"N={n:3d} dt={dt:6g} "
                  f"median CPU={statistics.median(measurements):.6f}s",
                  flush=True)

    csv_path = args.output / "baseline.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    def command_output(*command):
        result = subprocess.run(
            command, cwd=ROOT, text=True, capture_output=True)
        return result.stdout.strip() if result.returncode == 0 else None

    metadata = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": command_output("git", "rev-parse", "HEAD"),
        "git_describe": command_output("git", "describe",
                                       "--tags", "--always"),
        "git_status": command_output("git", "status", "--short"),
        "compiler": command_output("gfortran", "--version"),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "repetitions": args.repetitions,
        "meshes": args.meshes,
        "timesteps": args.timesteps,
        "duration": args.duration,
        "executable": str(executable),
        "timing": "Fortran CPU_TIME around INTEGRATE_TRANSIENT",
        "wall_timing": "Python perf_counter around subprocess",
        "numerical_baseline": "m9e-complete",
        "source": "existing INTEGRATE_TRANSIENT, unmodified",
    }
    (args.output / "baseline.json").write_text(
        json.dumps(metadata, indent=2) + "\n")
    print(f"Results: {csv_path}")

if __name__ == "__main__":
    main()
