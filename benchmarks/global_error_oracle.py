#!/usr/bin/env python3
"""M7c oracle experiment: separate defect estimation from error transport.

For y_num = G y_num_prev, y_exact = E y_exact_prev, and
error = y_num - y_exact, the exact identity is
 error_next = E error_prev + (G-E) y_num_prev.
The matrix exponential is an *oracle*, not a proposed production estimator.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import load_operators, ROOT

FIELDS = ['n', 'theta', 'dt_s', 'time_s', 'estimator',
          'actual_l2', 'estimated_l2', 'effectivity',
          'vector_relative_error', 'signed_alignment',
          'local_defect_relative_error', 'oracle_identity_relative_error']

def metrics(n, theta, h, time, name, actual, estimate, defect_rel, oracle_rel):
    a = np.linalg.norm(actual)
    b = np.linalg.norm(estimate)
    return dict(n=n, theta=theta, dt_s=h, time_s=time, estimator=name,
                actual_l2=a, estimated_l2=b,
                effectivity=b/a if a > 1e-14 else np.nan,
                vector_relative_error=np.linalg.norm(estimate-actual)/a
                if a > 1e-14 else np.nan,
                signed_alignment=np.dot(estimate,actual)/(a*b)
                if a*b > 1e-28 else np.nan,
                local_defect_relative_error=defect_rel,
                oracle_identity_relative_error=oracle_rel)

def run_case(j, n, theta, h, duration):
    ident = np.eye(j.shape[0]); steps = round(duration/h)
    if abs(steps*h-duration) > 1e-10:
        raise ValueError('duration must be divisible by timestep')
    lu = lu_factor(ident-theta*h*j)
    g = lu_solve(lu, ident+(1-theta)*h*j)
    e = expm(h*j)
    j2 = j @ j
    y = np.zeros(j.shape[0]); y[2*(n//2)+1] = 1.
    exact = y.copy()
    estimates = {name: np.zeros_like(y) for name in
                 ('oracle-exact', 'curvature-exact', 'resolvent-exact',
                  'curvature-theta', 'resolvent-theta')}
    rows = []
    for k in range(1,steps+1):
        defect = (g-e) @ y
        curvature = (theta-.5)*h*h*(j2 @ y)
        resolvent = lu_solve(lu, curvature)
        # The exact oracle identity uses E to propagate numerical-minus-exact
        # error, because its local defect is evaluated at the numerical state.
        for name, local in [('oracle-exact',defect),
                            ('curvature-exact',curvature),
                            ('resolvent-exact',resolvent)]:
            estimates[name] = e @ estimates[name] + local
        # These are intentionally heuristic, included to isolate the effect
        # of replacing E with G in the error propagation.
        for name, local in [('curvature-theta',curvature),
                            ('resolvent-theta',resolvent)]:
            estimates[name] = g @ estimates[name] + local
        y = g @ y; exact = e @ exact
        actual = y-exact
        oracle_rel = np.linalg.norm(estimates['oracle-exact']-actual) / max(
            np.linalg.norm(actual),1e-14)
        for name, estimate in estimates.items():
            local = defect if name == 'oracle-exact' else (
                curvature if 'curvature' in name else resolvent)
            defect_rel = np.linalg.norm(local-defect)/max(np.linalg.norm(defect),1e-14)
            rows.append(metrics(n,theta,h,k*h,name,actual,estimate,
                                defect_rel,oracle_rel))
    return rows

def write_csv(path, rows, fields):
    with path.open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=fields,lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'results-m7c')
    parser.add_argument('--duration',type=float,default=600.)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    all_rows=[]
    for n,j in sorted(load_operators().items()):
        for theta in (0.65,1.):
            for h in (7.5,15.,30.,60.):
                all_rows.extend(run_case(j,n,theta,h,args.duration))
    write_csv(args.output/'oracle_history.csv',all_rows,FIELDS)
    summary=[r for r in all_rows if abs(r['time_s']-args.duration)<1e-9]
    write_csv(args.output/'oracle_summary.csv',summary,FIELDS)
    print('n theta dt estimator effectivity vector-relative-error')
    for r in summary:
        print(f"{r['n']:2d} {r['theta']:.2f} {r['dt_s']:5.1f} "
              f"{r['estimator']:18s} {r['effectivity']:11.5g} "
              f"{r['vector_relative_error']:11.5g}")
    oracle_max=max(r['oracle_identity_relative_error'] for r in summary)
    print(f'Maximum oracle identity relative error: {oracle_max:.3e}')
    if not np.isfinite(oracle_max) or oracle_max > 1e-7:
        raise SystemExit('FAIL: exact-defect oracle does not close')
    print('PASS: exact-defect error-transport identity')
    print('Wrote',args.output/'oracle_summary.csv')

if __name__=='__main__':
    main()
