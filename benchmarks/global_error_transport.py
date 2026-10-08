#!/usr/bin/env python3
"""M7c exploratory Jacobian-based global error transport (linearized pipe).

Uses the M7b production-Jacobian export. The exact comparator is exp(TJ)u0,
not a finer numerical trajectory. No step doubling is performed.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve

ROOT = Path(__file__).resolve().parent

def load_operators():
    result = {}
    with (ROOT / 'results-m7b/modal_matrix.csv').open(newline='') as f:
        for row in csv.DictReader(f):
            n = int(row['n'])
            if n not in result:
                result[n] = np.zeros((2*n, 2*n), dtype=float)
            result[n][int(row['row'])-1, int(row['col'])-1] = float(row['operator'])
    return result

def experiment(j, h, theta, end_time, initial):
    steps = round(end_time / h)
    assert abs(steps*h-end_time) < 1e-10
    ident = np.eye(len(initial))
    lu = lu_factor(ident-theta*h*j)
    g = lu_solve(lu, ident+(1-theta)*h*j)
    # Exact semidiscrete reference; no temporal discretization error.
    exact_step = expm(h*j)
    y = initial.copy()
    reference = initial.copy()
    estimate = np.zeros_like(y)
    j2 = j @ j
    results = []
    for k in range(1, steps+1):
        # Heuristic resolvent-filtered leading defect; not an exact identity.
        tau = (theta-0.5)*h*h*lu_solve(lu, j2 @ y)
        estimate = g @ estimate + tau
        y = g @ y
        reference = exact_step @ reference
        actual = y-reference
        actual_norm = np.linalg.norm(actual)
        est_norm = np.linalg.norm(estimate)
        effectivity = est_norm/actual_norm if actual_norm > 1e-20 else float('nan')
        results.append(dict(time_s=k*h,actual_l2=actual_norm,
                            estimated_l2=est_norm,effectivity=effectivity,
                            actual_max=np.max(np.abs(actual)),
                            estimated_max=np.max(np.abs(estimate)),
                            vector_error_l2=np.linalg.norm(estimate-actual),
                            signed_alignment=(np.dot(estimate, actual)/
                              (est_norm*actual_norm) if est_norm*actual_norm>1e-30
                              else float('nan'))))
    return results

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'results-m7c')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    operators=load_operators()
    rows=[]
    for n,j in sorted(operators.items()):
        # A localized perturbation of a pressure-state coordinate. The
        # production operator's state ordering alternates flow and pressure.
        u0=np.zeros(2*n)
        u0[2*(n//2)+1]=1.0
        for theta in (0.65, 1.0):
            for h in (7.5, 15., 30., 60.):
                history=experiment(j,h,theta,600.,u0)
                for row in history:
                    rows.append(dict(n=n,dt_s=h,theta=theta,**row))
    with (args.output/'transport_history.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    summary=[]
    for n in sorted(operators):
        for theta in (0.65,1.0):
            for h in (7.5,15.,30.,60.):
                r=[r for r in rows if r['n']==n and r['theta']==theta and r['dt_s']==h]
                final=r[-1]
                peak=max(r,key=lambda x:x['actual_l2'])
                summary.append(dict(n=n,dt_s=h,theta=theta,
                    final_effectivity=final['effectivity'],
                    final_vector_relative_error=(final['vector_error_l2']/final['actual_l2']
                       if final['actual_l2']>1e-20 else float('nan')),
                    peak_error_time_s=peak['time_s'],
                    peak_error_effectivity=peak['effectivity'],
                    max_actual_l2=peak['actual_l2'],
                    max_estimated_l2=max(x['estimated_l2'] for x in r)))
    with (args.output/'transport_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader();writer.writerows(summary)
    print('n theta dt final-effectivity peak-effectivity final-vector-rel-error')
    for x in summary:
        print(f"{x['n']:2d} {x['theta']:.2f} {x['dt_s']:5.1f} "
              f"{x['final_effectivity']:10.4f} {x['peak_error_effectivity']:10.4f} "
              f"{x['final_vector_relative_error']:10.4f}")
    print('Wrote',args.output)

if __name__=='__main__':
    main()
