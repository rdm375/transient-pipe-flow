#!/usr/bin/env python3
"""Numerical regression checks for the standalone C++ simulator."""

import argparse
import math
import subprocess
import sys


REFERENCES = {
    80: {
        "final_outlet_pressure_pa": 7846825.3670041412,
        "final_inlet_flow_kg_s": 103.19772530073323,
        "final_outlet_flow_kg_s": 105.0,
        "final_linepack_kg": 4799630.9833060475,
        "max_mass_defect_kg_s": 3.9544278962466706e-10,
        "total_newton_iterations": 160,
    },
    8000: {
        "final_outlet_pressure_pa": 7843657.8677442567,
        "final_inlet_flow_kg_s": 104.99999993783844,
        "final_outlet_flow_kg_s": 105.0,
        "final_linepack_kg": 4798420.4378263094,
        "max_mass_defect_kg_s": 6.2161561231732776e-08,
        "total_newton_iterations": 1049,
    },
}


def parse_output(output):
    values = {}

    for line in output.splitlines():
        if "=" not in line:
            continue

        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()

    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("steps", type=int)
    args = parser.parse_args()

    if args.steps not in REFERENCES:
        parser.error("No reference for this step count")

    result = subprocess.run(
        [args.executable, "100", str(args.steps)],
        capture_output=True,
        text=True,
        check=True,
    )

    values = parse_output(result.stdout)
    reference = REFERENCES[args.steps]

    failures = []

    for key, expected in reference.items():
        if key not in values:
            failures.append(f"Missing output: {key}")
            continue

        if key == "total_newton_iterations":
            actual = int(values[key])

            if actual != expected:
                failures.append(
                    f"{key}: actual={actual}, expected={expected}"
                )

            continue

        actual = float(values[key])

        if not math.isfinite(actual):
            failures.append(f"{key}: nonfinite value {actual}")
            continue

        if key == "max_mass_defect_kg_s":
            # The defect is a small difference of large quantities.
            # Use an absolute bound, not a relative tolerance.
            tolerance = 1.0e-9
        else:
            tolerance = 5.0e-13 * abs(expected)

        error = abs(actual - expected)

        if error > tolerance:
            failures.append(
                f"{key}: actual={actual:.17g}, "
                f"expected={expected:.17g}, "
                f"error={error:.6e}, "
                f"tolerance={tolerance:.6e}"
            )

    if failures:
        for failure in failures:
            print("FAIL:", failure, file=sys.stderr)
        return 1

    print(
        f"PASS C++ numerical regression: "
        f"cells=100 steps={args.steps}"
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
