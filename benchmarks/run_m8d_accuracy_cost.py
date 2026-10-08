#!/usr/bin/env python3
"""M8D: fixed versus adaptive accuracy/cost comparison."""
import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path

import numpy as np

from run_adaptive_inventory_m8 import SOURCES, read


def execute(command, repeats):
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.DEVNULL,
        )
        times.append(time.perf_counter() - start)
    return statistics.median(times)


def errors(path, reference):
    result = read(path)
    t = result["time_s"]

    mass_reference = np.interp(
        t, reference["time_s"], reference["inventory_kg"]
    )
    flow_reference = np.interp(
        t, reference["time_s"], reference["imbalance_kg_s"]
    )

    mass_error = result["inventory_kg"] - mass_reference
    flow_error = result["imbalance_kg_s"] - flow_reference

    return {
        "steps": len(t) - 1,
        "inventory_error_kg": float(
            np.max(np.abs(mass_error - mass_error[0]))
        ),
        "imbalance_error_kg_s": float(np.max(np.abs(flow_error))),
        "conservation_defect_kg": float(
            np.max(np.abs(result["step_balance_defect_kg"][1:]))
        ),
        "accepted_newton_iterations": int(
            np.sum(result["newton_iterations"])
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, nargs="+", default=[20, 40])
    parser.add_argument(
        "--fixed-dt", type=float, nargs="+",
        default=[7.5, 15, 30, 60, 120],
    )
    parser.add_argument(
        "--tol", type=float, nargs="+",
        default=[0.003, 0.01, 0.03],
    )
    parser.add_argument("--duration", type=float, default=3600)
    parser.add_argument("--reference-dt", type=float, default=3.75)
    parser.add_argument("--hmax", type=float, default=120)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument(
        "--outdir", type=Path,
        default=Path("benchmarks/results-m8/m8d-accuracy-cost"),
    )
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be positive")

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    fixed_exe = outdir / "fixed_solver"
    adaptive_exe = outdir / "adaptive_solver"

    flags = ["gfortran", "-O2", "-ffixed-line-length-none"]

    subprocess.run(
        [
            *flags, "-o", str(fixed_exe), *SOURCES,
            "benchmarks/inventory_history_m7c18a.f",
        ],
        check=True,
    )
    subprocess.run(
        [
            *flags, "-o", str(adaptive_exe), *SOURCES,
            "benchmarks/adaptive_inventory_m8.f",
        ],
        check=True,
    )

    rows = []

    for n in args.n:
        reference_path = outdir / f"reference_n{n}.csv"

        subprocess.run(
            [
                str(fixed_exe),
                str(args.reference_dt),
                str(args.duration),
                str(reference_path),
                str(n),
            ],
            check=True,
        )

        reference = read(reference_path)

        for method, parameters in (
            ("fixed", args.fixed_dt),
            ("adaptive", args.tol),
        ):
            for parameter in parameters:
                path = outdir / f"{method}_n{n}_{parameter:g}.csv"

                if method == "fixed":
                    command = [
                        str(fixed_exe),
                        str(parameter),
                        str(args.duration),
                        str(path),
                        str(n),
                    ]
                else:
                    command = [
                        str(adaptive_exe),
                        str(path),
                        str(n),
                        str(args.duration),
                        str(parameter),
                        str(args.hmax),
                    ]

                elapsed = execute(command, args.repeats)
                metrics = errors(path, reference)

                rejected = 0
                if method == "adaptive":
                    data = read(path)
                    rejected = int(
                        np.sum(data["rejects_before_accept"])
                    )

                row = {
                    "n": n,
                    "method": method,
                    "parameter": parameter,
                    "median_wall_s": elapsed,
                    "accepted_steps": metrics["steps"],
                    "rejected_steps": rejected,
                    "attempts": metrics["steps"] + rejected,
                    **{
                        k: v for k, v in metrics.items()
                        if k != "steps"
                    },
                }

                rows.append(row)

                print(
                    f"n={n:3d} {method:8s} "
                    f"parameter={parameter:8g} "
                    f"time={elapsed:.6f}s "
                    f"inventory={metrics['inventory_error_kg']:.5g}kg "
                    f"imbalance={metrics['imbalance_error_kg_s']:.5g}kg/s",
                    flush=True,
                )

    summary = outdir / "summary.csv"
    with summary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("\nAccuracy-matched comparisons")
    print("----------------------------")

    for n in args.n:
        fixed = [
            r for r in rows
            if r["n"] == n and r["method"] == "fixed"
        ]
        adaptive = [
            r for r in rows
            if r["n"] == n and r["method"] == "adaptive"
        ]

        for a in adaptive:
            eligible = [
                f for f in fixed
                if f["inventory_error_kg"]
                   <= a["inventory_error_kg"]
                and f["imbalance_error_kg_s"]
                   <= a["imbalance_error_kg_s"]
            ]

            if not eligible:
                print(
                    f"n={n} tol={a['parameter']:g}: "
                    "no fixed-step run satisfies both error limits"
                )
                continue

            best = min(eligible, key=lambda r: r["median_wall_s"])
            speedup = best["median_wall_s"] / a["median_wall_s"]

            print(
                f"n={n} tol={a['parameter']:g} "
                f"fixed_dt={best['parameter']:g} "
                f"speedup={speedup:.3f}x"
            )

    print(f"\nWrote {summary}")
    print(
        "Timing includes process startup and CSV output. "
        "Newton counts exclude rejected attempts."
    )
    print(
        "Accuracy is measured against a finite-step reference "
        "and is not certified."
    )


if __name__ == "__main__":
    main()
