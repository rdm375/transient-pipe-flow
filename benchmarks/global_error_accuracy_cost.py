#!/usr/bin/env python3
"""M7c-5: residual-driven CN versus TR-BDF2 error transport.

The matrix exponential is used ONLY to calculate the reference error.
Timing is diagnostic, not a production-solver overhead estimate.
"""
import argparse
import csv
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import ROOT, load_operators

GAMMA = 2.0 - np.sqrt(2.0)
FIELDS = ('n', 'production_theta', 'dt_s', 'method', 'substeps',
          'effectivity', 'vector_relative_error', 'signed_alignment',
          'factor_seconds', 'solve_seconds', 'transport_seconds',
          'factorizations', 'linear_solves', 'oracle_relative_error')


def run_case(j, n, theta, h, duration, counts):
    steps = round(duration / h)
    if abs(steps*h-duration) > 1e-9:
        raise ValueError('Duration must be divisible by dt')
    eye = np.eye(len(j))
    lu_prod = lu_factor(eye-theta*h*j)
    g = lu_solve(lu_prod, eye+(1-theta)*h*j)
    exact_step = expm(h*j)  # validation ONLY
    configs = {}
    for method in ('CN', 'TR-BDF2'):
        for m in counts:
            ds = h/m
            t0 = perf_counter()
            if method == 'CN':
                lu1 = lu_factor(eye-0.5*ds*j)
                data = (lu1, None, eye+0.5*ds*j)
                factors = 1
            else:
                # gamma = 2-sqrt(2): both stage matrices have identical
                # coefficient gamma/2 = (1-gamma)/(2-gamma).
                alpha = GAMMA/2
                lu1 = lu_factor(eye-alpha*ds*j)
                data = (lu1, None, eye+alpha*ds*j)
                factors = 1
            factor_seconds = perf_counter()-t0
            configs[(method,m)] = dict(data=data, err=np.zeros(len(j)),
                factor_seconds=factor_seconds, solve_seconds=0.,
                transport_seconds=0., linear_solves=0, factorizations=factors)

    y = np.zeros(len(j))
    y[2*(n//2)+1] = 1.
    reference = y.copy()
    oracle = np.zeros_like(y)
    for _ in range(steps):
        next_y = g@y
        v = (next_y-y)/h
        r0 = v-j@y
        dr = -j@(next_y-y)
        for (method,m), c in configs.items():
            ds = h/m
            lu, _, left = c['data']
            err = c['err']
            begin = perf_counter()
            for k in range(m):
                frac = k/m
                r_start = r0+frac*dr
                r_end = r0+(frac+1/m)*dr
                if method == 'CN':
                    rhs = left@err+0.5*ds*(r_start+r_end)
                    t0=perf_counter()
                    err = lu_solve(lu,rhs)
                    c['solve_seconds'] += perf_counter()-t0
                    c['linear_solves'] += 1
                else:
                    r_mid = r0+(frac+GAMMA/m)*dr
                    rhs1 = left@err + (GAMMA*ds/2)*(r_start+r_mid)
                    t0=perf_counter()
                    middle = lu_solve(lu,rhs1)
                    c['solve_seconds'] += perf_counter()-t0
                    c['linear_solves'] += 1
                    # Variable-step BDF2 from 0, gamma*ds, ds.
                    # (I - alpha*ds*J) err_end = weighted states + alpha*ds*r_end
                    rhs2 = (middle/(GAMMA*(2-GAMMA))
                            - ((1-GAMMA)**2)/(GAMMA*(2-GAMMA))*err
                            + ((1-GAMMA)/(2-GAMMA))*ds*r_end)
                    t0=perf_counter()
                    err = lu_solve(lu,rhs2)
                    c['solve_seconds'] += perf_counter()-t0
                    c['linear_solves'] += 1
            c['transport_seconds'] += perf_counter()-begin
            c['err'] = err
        oracle = exact_step@oracle + next_y-exact_step@y
        y = next_y
        reference = exact_step@reference
    actual = y-reference
    anorm = np.linalg.norm(actual)
    oracle_rel = np.linalg.norm(oracle-actual)/max(anorm,1e-14)
    rows=[]
    for (method,m), c in configs.items():
        err=c['err']; enorm=np.linalg.norm(err)
        rows.append(dict(n=n, production_theta=theta, dt_s=h,
            method=method, substeps=m, effectivity=enorm/anorm,
            vector_relative_error=np.linalg.norm(err-actual)/anorm,
            signed_alignment=np.dot(err,actual)/(enorm*anorm) if enorm*anorm>1e-28 else np.nan,
            factor_seconds=c['factor_seconds'],solve_seconds=c['solve_seconds'],
            transport_seconds=c['transport_seconds'],
            factorizations=c['factorizations'],linear_solves=c['linear_solves'],
            oracle_relative_error=oracle_rel))
    return rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--duration',type=float,default=600.)
    parser.add_argument('--substeps',type=int,nargs='+',default=[8,16,32,64,128,256])
    parser.add_argument('--grids',type=int,nargs='+',default=[20,40])
    args=parser.parse_args()
    if any(m<1 for m in args.substeps):
        parser.error('substeps must be positive')
    rows=[]
    for n,j in sorted(load_operators().items()):
        if n not in args.grids: continue
        for theta in (0.65,1.):
            for h in (30.,60.):
                rows.extend(run_case(j,n,theta,h,args.duration,args.substeps))
    if not rows: parser.error('no matching operators')
    out=ROOT/'results-m7c';out.mkdir(parents=True,exist_ok=True)
    path=out/'accuracy_cost_summary.csv'
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(rows)
    max_oracle=max(r['oracle_relative_error'] for r in rows)
    print('Max oracle relative error: %.3e'%max_oracle)
    if not np.isfinite(max_oracle) or max_oracle>1e-9:
        raise SystemExit('FAIL: oracle transport identity')
    print('PASS: oracle transport identity')
    print('n prod-theta dt method substeps effectivity vector-error solves solve-ms')
    for r in rows:
        print('%2d %.2f %4.0f %-7s %3d %11.6g %12.6g %5d %9.2f'%(
            r['n'],r['production_theta'],r['dt_s'],r['method'],r['substeps'],
            r['effectivity'],r['vector_relative_error'],r['linear_solves'],
            1e3*r['solve_seconds']))
    print('Wrote',path)

if __name__=='__main__': main()
