#!/usr/bin/env python3
"""M7c-13: compare augmented residual action with reduced J-only recurrence.

Error convention: numerical minus exact.
Direct identity: e_next = exp(hJ) (e_current - y_current) + y_next.
No step doubling. Dense expm is used only for independent validation.
Fixed Arnoldi comparisons use identical Krylov dimensions (not equal accuracy).
"""
import argparse
from pathlib import Path
import csv
import time
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from scipy.sparse import bmat, csc_matrix, eye
from scipy.sparse.linalg import expm_multiply
from global_error_transport import ROOT, load_operators
from global_error_adaptive_krylov import initial_state
from global_error_budgeting import candidates

def relative(a, b):
    return np.linalg.norm(a-b) / max(np.linalg.norm(b), 1e-14)

def run(j, n, theta, h, duration, probe, seed, dims):
    size = len(j)
    ident = np.eye(size)
    g = lu_solve(lu_factor(ident-theta*h*j),
                 ident+(1-theta)*h*j)
    e = expm(h*j)  # reference only
    js = csc_matrix(j)
    z = csc_matrix((size, size))
    a = bmat([[js, eye(size,format='csc'), z],
              [z, z, eye(size,format='csc')],
              [z, z, z]], format='csc')
    steps = round(duration/h)
    if not np.isclose(steps*h, duration):
        raise ValueError("duration must be divisible by dt")
    y0 = initial_state(j,n,probe,seed)
    y = y0.copy()
    exact = y0.copy()
    ref = np.zeros(size)
    aug_high = np.zeros(size)
    direct_high = np.zeros(size)
    aug_fixed = {m:np.zeros(size) for m in dims}
    direct_fixed = {m:np.zeros(size) for m in dims}
    elapsed = {'aug_high':0., 'direct_high':0.}
    for m in dims:
        elapsed[f'aug_{m}']=0.
        elapsed[f'direct_{m}']=0.
    step_rows=[]
    for k in range(steps):
        ynext = g@y
        r0 = (ynext-y)/h-j@y
        r1 = (ynext-y)/h-j@ynext
        w = (r1-r0)/h
        ref = e@ref + (g-e)@y
        exact = e@exact

        vaug=np.concatenate((aug_high,r0,w))
        t0=time.perf_counter()
        aug_high = expm_multiply(h*a, vaug)[:size]
        elapsed['aug_high']+=time.perf_counter()-t0

        # Reduced action only on J. The additive y_next is exact up to
        # floating point; cancellation can occur in the final addition.
        vdir=direct_high-y
        t0=time.perf_counter()
        direct_high=expm_multiply(h*js,vdir)+ynext
        elapsed['direct_high']+=time.perf_counter()-t0

        for m in dims:
            vaug=np.concatenate((aug_fixed[m],r0,w))
            t0=time.perf_counter()
            aug_fixed[m]=candidates(a,vaug,h,[m],size)[m]
            elapsed[f'aug_{m}']+=time.perf_counter()-t0
            vdir=direct_fixed[m]-y
            t0=time.perf_counter()
            direct_fixed[m]=candidates(js,vdir,h,[m],size)[m]+ynext
            elapsed[f'direct_{m}']+=time.perf_counter()-t0

        # Estimate sensitivity of direct reconstruction to cancellation:
        # ratio of sum of term magnitudes to magnitude of their sum.
        propagated=expm_multiply(h*js, direct_high-y)
        cancellation=(np.linalg.norm(propagated)+np.linalg.norm(ynext)) / max(
            np.linalg.norm(propagated+ynext),1e-30)
        step_rows.append(dict(n=n,theta=theta,dt_s=h,probe=probe,step=k+1,
                              time_s=(k+1)*h,true_error_norm=np.linalg.norm(ref),
                              direct_high_rel=relative(direct_high,ref),
                              augmented_high_rel=relative(aug_high,ref),
                              direct_cancellation_ratio=cancellation))
        y=ynext
    true=y-exact
    denom=max(np.linalg.norm(true),1e-14)
    rows=[]
    for method, estimate, key, jv in [
        ('augmented_scipy',aug_high,'aug_high',None),
        ('reduced_scipy',direct_high,'direct_high',None),
    ]+[(f'augmented_arnoldi_m{m}',aug_fixed[m],f'aug_{m}',steps*m)
       for m in dims]+[(f'reduced_arnoldi_m{m}',direct_fixed[m],f'direct_{m}',steps*m)
       for m in dims]:
        rows.append(dict(n=n,theta=theta,dt_s=h,probe=probe,method=method,
                         steps=steps,jacobian_products=jv if jv is not None else '',
                         elapsed_ms=1000*elapsed[key],
                         vector_relative_error=np.linalg.norm(estimate-true)/denom,
                         absolute_error=np.linalg.norm(estimate-true),
                         true_error_norm=np.linalg.norm(true),
                         oracle_recurrence_error=np.linalg.norm(ref-true)/denom,
                         estimate_norm=np.linalg.norm(estimate)))
    return rows,step_rows

def write_csv(path, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print('Wrote',path)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,nargs='+',default=[40])
    p.add_argument('--theta',type=float,nargs='+',default=[.65])
    p.add_argument('--dt',type=float,nargs='+',default=[60.])
    p.add_argument('--probes',nargs='+',default=['fast','slow','impulse'],
                   choices=['fast','slow','impulse','random'])
    p.add_argument('--dimensions',type=int,nargs='+',default=[8,12,16,20,24,28,32])
    p.add_argument('--duration',type=float,default=600.)
    p.add_argument('--seed',type=int,default=1729)
    p.add_argument('--output-dir',type=Path,default=ROOT/'results-m7c')
    args=p.parse_args()
    if any(m<1 for m in args.dimensions): p.error('dimensions must be positive')
    rows=[]; steps=[]
    operators=load_operators()
    for n in args.n:
        for theta in args.theta:
            for dt in args.dt:
                for probe in args.probes:
                    a,b=run(operators[n],n,theta,dt,args.duration,probe,args.seed,
                            sorted(set(args.dimensions)))
                    rows.extend(a); steps.extend(b)
                    for row in a:
                        print(f'{n} theta={theta:g} dt={dt:g} {probe:7s} '
                              f'{row["method"]:24s} '
                              f'err={row["vector_relative_error"]:.5g} '
                              f'Jv={str(row["jacobian_products"]):>4s} '
                              f'ms={row["elapsed_ms"]:.1f}',flush=True)
    write_csv(args.output_dir/'reduced_recurrence_comparison.csv',rows)
    write_csv(args.output_dir/'reduced_recurrence_steps.csv',steps)
    print('Fixed-Arnoldi Jv counts include one J matvec per basis dimension;')
    print('SciPy Jv counts are not instrumented. Timings include Python overhead.')
    print('The direct formula can suffer cancellation: inspect step diagnostics.')
if __name__=='__main__':
    main()
