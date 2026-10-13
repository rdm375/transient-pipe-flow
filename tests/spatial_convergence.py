#!/usr/bin/env python3
"""Spatial-convergence regression for the native C++ transient solver."""

import argparse
import math
import subprocess


REFERENCE = {
    100: {
        "final_outlet_pressure_pa": 7846825.3670041412,
        "final_inlet_flow_kg_s": 103.19772530073323,
        "final_linepack_kg": 4799630.9833060475,
    },
    200: {
        "final_outlet_pressure_pa": 7846825.2498928932,
        "final_inlet_flow_kg_s": 103.19773419864673,
        "final_linepack_kg": 4799630.9687672555,
    },
    400: {
        "final_outlet_pressure_pa": 7846825.2206157222,
        "final_inlet_flow_kg_s": 103.19773642323406,
        "final_linepack_kg": 4799630.9651325513,
    },
}

FIELDS = (
    "final_outlet_pressure_pa",
    "final_inlet_flow_kg_s",
    "final_linepack_kg",
)


def run_simulator(executable, cells):
    completed = subprocess.run(
        [executable, str(cells), "80"],
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )

    values = {}

    for line in completed.stdout.splitlines():
        if "=" not in line:
            continue

        key, value = line.split("=", 1)

        if key == "backend":
            values[key] = value
        else:
            values[key] = float(value)

    required = (
        "cells",
        "steps",
        "final_outlet_pressure_pa",
        "final_inlet_flow_kg_s",
        "final_linepack_kg",
        "max_mass_defect_kg_s",
        "total_newton_iterations",
    )

    for key in required:
        if key not in values:
            raise AssertionError(
                f"n={cells}: missing output field {key}"
            )

    assert values["backend"] == "cpp20"
    assert values["cells"] == cells
    assert values["steps"] == 80

    for key in required:
        if not math.isfinite(values[key]):
            raise AssertionError(
                f"n={cells}: nonfinite {key}={values[key]}"
            )

    return values


def convergence_order(u100, u200, u400):
    coarse = abs(u100 - u200)
    fine = abs(u200 - u400)

    if fine == 0.0:
        raise AssertionError(
            "Fine-grid difference is zero; cannot estimate order"
        )

    if coarse == 0.0:
        raise AssertionError(
            "Coarse-grid difference is zero; cannot estimate order"
        )

    return math.log2(coarse / fine)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "executable",
        help="Path to the pipe_sim_cpp executable",
    )

    args = parser.parse_args()

    results = {
        n: run_simulator(args.executable, n)
        for n in (100, 200, 400)
    }

    print("Spatial convergence: 100 -> 200 -> 400 cells")

    for field in FIELDS:
        values = [
            results[n][field]
            for n in (100, 200, 400)
        ]

        order = convergence_order(*values)

        print(f"{field}:")
        print(f"  n=100: {values[0]:.17g}")
        print(f"  n=200: {values[1]:.17g}")
        print(f"  n=400: {values[2]:.17g}")
        print(f"  observed_order={order:.8f}")

        # The current model exhibits approximately second-order
        # spatial convergence for these three output quantities.
        if not 1.8 <= order <= 2.2:
            raise AssertionError(
                f"{field}: unexpected convergence order {order}"
            )

        # Guard against a common-mode shift that preserves the
        # observed convergence order.
        #
        # These absolute tolerances are deliberately wider than
        # floating-point roundoff, but small relative to the
        # physical scale of each quantity.
        atol = {
            "final_outlet_pressure_pa": 0.02,
            "final_inlet_flow_kg_s": 2.0e-6,
            "final_linepack_kg": 0.005,
        }[field]

        for n in (100, 200, 400):
            actual = results[n][field]
            expected = REFERENCE[n][field]

            if abs(actual - expected) > atol:
                raise AssertionError(
                    f"n={n}, {field}: "
                    f"actual={actual:.17g}, "
                    f"reference={expected:.17g}, "
                    f"tolerance={atol}"
                )

    for n, result in results.items():
        defect = result["max_mass_defect_kg_s"]

        print(
            f"n={n}: "
            f"mass_defect={defect:.6e} kg/s, "
            f"newton_iterations="
            f"{int(result['total_newton_iterations'])}"
        )

        if not 0.0 <= defect < 1.0e-8:
            raise AssertionError(
                f"n={n}: excessive mass conservation defect"
            )

        iterations = result["total_newton_iterations"]

        # The reference case uses 160 iterations for 80 steps.
        # Allow modest variation rather than requiring exact equality.
        if not 80 <= iterations <= 240:
            raise AssertionError(
                f"n={n}: unexpected Newton iteration count"
            )

    print("PASS spatial convergence regression")


if __name__ == "__main__":
    main()
