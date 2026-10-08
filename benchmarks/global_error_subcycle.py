#!/usr/bin/env python3
"""M7c-4: subcycled residual-driven global error transport."""
import csv
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import ROOT, load_operators

FIELDS = [
    "n", "theta", "dt_s", "transport_theta", "substeps",
    "effectivity", "vector_relative_error", "signed_alignment",
    "oracle_relative_error"
]

def experiment(j, n, theta, h):
    I = np.eye(len(j))
    G = np.linalg.solve(I-theta*h*j, I+(1-theta)*h*j)
    E = expm(h*j)  # Reference only; not used by estimators.

    y = np.zeros(len(j))
    y[2*(n//2)+1] = 1.0
    reference = y.copy()
    oracle = np.zeros_like(y)

    configurations = {}
    for a in (0.5, 0.65, 1.0):
        for m in (1, 2, 4, 8, 16, 32, 64):
            d = h/m
            configurations[a,m] = (
                lu_factor(I-a*d*j),
                I+(1-a)*d*j,
                d,
                np.zeros_like(y)
            )

    for _ in range(round(600/h)):
        ynext = G@y
        v = (ynext-y)/h
        r0 = v-j@y
        dr = -j@(ynext-y)

        for (a,m), (lu,L,d,err) in configurations.items():
            for k in range(m):
                ra = r0+(k/m)*dr
                rb = r0+((k+1)/m)*dr
                err = lu_solve(
                    lu, L@err+d*((1-a)*ra+a*rb)
                )
            configurations[a,m] = (lu,L,d,err)

        oracle = E@oracle+ynext-E@y
        y = ynext
        reference = E@reference

    actual = y-reference
    norm = np.linalg.norm(actual)
    oracle_error = np.linalg.norm(oracle-actual)/max(norm,1e-14)
    rows = []

    for (a,m), (_,_,_,err) in configurations.items():
        estimated = np.linalg.norm(err)
        rows.append(dict(
            n=n, theta=theta, dt_s=h,
            transport_theta=a, substeps=m,
            effectivity=estimated/norm,
            vector_relative_error=np.linalg.norm(err-actual)/norm,
            signed_alignment=(
                np.dot(err,actual)/(estimated*norm)
                if estimated*norm>1e-28 else np.nan
            ),
            oracle_relative_error=oracle_error
        ))
    return rows

def main():
    rows = []
    for n,j in sorted(load_operators().items()):
        for theta in (0.65,1.0):
            for h in (7.5,15.0,30.0,60.0):
                rows.extend(experiment(j,n,theta,h))

    output = ROOT/"results-m7c"/"subcycle_summary.csv"
    output.parent.mkdir(parents=True,exist_ok=True)

    with output.open("w",newline="") as f:
        writer = csv.DictWriter(f,FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    oracle_max = max(r["oracle_relative_error"] for r in rows)
    print(f"Maximum oracle relative error: {oracle_max:.3e}")
    assert oracle_max < 1e-9
    print("PASS: oracle transport identity")

    print("n theta dt transport substeps effectivity vector-error")
    for r in rows:
        if r["dt_s"] in (30.0,60.0) and r["substeps"] in (1,8,64):
            print(
                f'{r["n"]:2d} {r["theta"]:.2f} '
                f'{r["dt_s"]:4.0f} {r["transport_theta"]:.2f} '
                f'{r["substeps"]:3d} {r["effectivity"]:10.5g} '
                f'{r["vector_relative_error"]:10.5g}'
            )
    print("Wrote",output)

if __name__ == "__main__":
    main()
