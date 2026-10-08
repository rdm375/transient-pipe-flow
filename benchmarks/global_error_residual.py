#!/usr/bin/env python3
"""M7c: residual reconstruction and global error transport."""
import csv
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import ROOT, load_operators

ORDERS = (2, 4, 8, 16)
FIELDS = (
    "n", "theta", "dt_s", "time_s", "estimator",
    "quadrature_points", "actual_l2", "estimated_l2",
    "effectivity", "vector_relative_error", "signed_alignment",
    "local_defect_relative_error", "oracle_identity_relative_error",
    "residual_identity_relative_error",
)


def residual_map(j, g, h, q):
    """Integrate exp((h-s)J) times the linear reconstruction residual."""
    ident = np.eye(len(j))
    delta = g - ident
    nodes, weights = leggauss(q)
    result = np.zeros_like(j)

    for x, w in zip(nodes, weights):
        s = h * (x + 1) / 2
        residual = delta / h - j @ (ident + (s / h) * delta)
        result += (h * w / 2) * (expm((h - s) * j) @ residual)

    return result


def run_case(j, n, theta, h, duration=600.0):
    steps = round(duration / h)
    assert abs(steps * h - duration) < 1e-10

    ident = np.eye(len(j))
    lu = lu_factor(ident - theta * h * j)
    g = lu_solve(lu, ident + (1 - theta) * h * j)
    e = expm(h * j)
    defect_map = g - e

    maps = {q: residual_map(j, g, h, q) for q in ORDERS}
    map_error = np.linalg.norm(maps[16] - defect_map) / max(
        np.linalg.norm(defect_map), 1e-14
    )

    y = np.zeros(len(j))
    y[2 * (n // 2) + 1] = 1.0
    reference = y.copy()

    estimates = {"oracle-exact": np.zeros_like(y)}
    for q in ORDERS:
        estimates[f"residual-q{q}-exact"] = np.zeros_like(y)
        estimates[f"residual-q{q}-theta"] = np.zeros_like(y)

    rows = []
    max_oracle_error = 0.0

    for k in range(1, steps + 1):
        prior = y.copy()
        defect = defect_map @ prior

        estimates["oracle-exact"] = (
            e @ estimates["oracle-exact"] + defect
        )

        for q in ORDERS:
            local = maps[q] @ prior
            exact_name = f"residual-q{q}-exact"
            theta_name = f"residual-q{q}-theta"

            estimates[exact_name] = (
                e @ estimates[exact_name] + local
            )
            estimates[theta_name] = (
                g @ estimates[theta_name] + local
            )

        y = g @ y
        reference = e @ reference
        actual = y - reference
        a = np.linalg.norm(actual)

        oracle_error = (
            np.linalg.norm(estimates["oracle-exact"] - actual)
            / max(a, 1e-14)
        )
        max_oracle_error = max(max_oracle_error, oracle_error)

        for name, estimated in estimates.items():
            q = 0 if name == "oracle-exact" else int(
                name.split("-")[1][1:]
            )
            local = defect if q == 0 else maps[q] @ prior
            b = np.linalg.norm(estimated)

            rows.append(dict(
                n=n, theta=theta, dt_s=h, time_s=k*h,
                estimator=name, quadrature_points=q,
                actual_l2=a, estimated_l2=b,
                effectivity=b/a if a > 1e-14 else np.nan,
                vector_relative_error=(
                    np.linalg.norm(estimated-actual)/a
                    if a > 1e-14 else np.nan
                ),
                signed_alignment=(
                    np.dot(estimated,actual)/(a*b)
                    if a*b > 1e-28 else np.nan
                ),
                local_defect_relative_error=(
                    np.linalg.norm(local-defect)
                    / max(np.linalg.norm(defect),1e-14)
                ),
                oracle_identity_relative_error=oracle_error,
                residual_identity_relative_error=map_error,
            ))

    return rows, max_oracle_error, map_error


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=FIELDS, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    output = ROOT / "results-m7c"
    output.mkdir(parents=True, exist_ok=True)

    rows = []
    oracle_errors = []
    residual_errors = []

    for n, j in sorted(load_operators().items()):
        for theta in (0.65, 1.0):
            for h in (7.5, 15.0, 30.0, 60.0):
                result, oracle_error, residual_error = run_case(
                    j, n, theta, h
                )
                rows.extend(result)
                oracle_errors.append(oracle_error)
                residual_errors.append(residual_error)

    summary = [
        row for row in rows
        if abs(row["time_s"] - 600.0) < 1e-9
    ]

    write_csv(output / "residual_history.csv", rows)
    write_csv(output / "residual_summary.csv", summary)

    print("n theta dt estimator effectivity vector-relative-error")
    for r in summary:
        print(
            f'{r["n"]:2d} {r["theta"]:.2f} {r["dt_s"]:5.1f} '
            f'{r["estimator"]:20s} {r["effectivity"]:11.5g} '
            f'{r["vector_relative_error"]:11.5g}'
        )

    oracle_max = max(oracle_errors)
    residual_max = max(residual_errors)

    print(f"Max oracle closure relative error: {oracle_max:.3e}")
    print(
        "Max highest-order residual operator relative error: "
        f"{residual_max:.3e}"
    )

    if oracle_max > 1e-9 or residual_max > 1e-8:
        raise SystemExit("FAIL: residual/oracle identity")

    print("PASS: oracle and residual identities")


if __name__ == "__main__":
    main()
