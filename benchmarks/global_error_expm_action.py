#!/usr/bin/env python3
"""M7c-6b: independently verify residual-driven error transport by expm_multiply.

The candidate uses an augmented sparse exponential action, never a dense
matrix exponential. Dense expm is used ONLY to construct the reference.
Error convention: numerical minus exact.
"""
import argparse
import csv
import time
from pathlib import Path
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from scipy.sparse import csc_matrix, bmat
from scipy.sparse.linalg import expm_multiply, LinearOperator
from global_error_transport import ROOT, load_operators

FIELDS = ('n', 'production_theta', 'dt_s', 'steps', 'method', 'effectivity',
          'vector_relative_error', 'oracle_relative_error', 'elapsed_ms',
          'jacobian_products', 'relative_local_defect_error')


def run_case(j, n, theta, h, duration):
    size = j.shape[0]
    ident = np.eye(size)
    lu = lu_factor(ident - theta * h * j)
    g = lu_solve(lu, ident + (1.0-theta)*h*j)
    e = expm(h*j)  # validation only
    # The autonomous augmented equation is
    # [z'; a'; b'] = [J z+a; b; 0].
    # a(0)=r0 and b(0)=(r1-r0)/h, so a(t)=r(t).
    js = csc_matrix(j)
    # Instead use the fixed 3-block system [z; v; w], where
    # z' = Jz + v, v'=w, w'=0 and v,w are size-dimensional vectors.
    # This is a 3*size operator but its sparse structure is simple.
    zeros = csc_matrix((size, size))
    eye = csc_matrix(ident)
    augmented = bmat([[js, eye, zeros],
                      [zeros, zeros, eye],
                      [zeros, zeros, zeros]], format='csc')
    counts = {"columns": 0}
    def mv(x):
        counts["columns"] += 1
        return augmented @ x
    def mm(x):
        counts["columns"] += x.shape[1]
        return augmented @ x
    def rmv(x):
        return augmented.T @ x

    def rmm(x):
        return augmented.T @ x

    action = LinearOperator(
        augmented.shape,
        matvec=mv,
        matmat=mm,
        rmatvec=rmv,
        rmatmat=rmm,
        dtype=np.float64,
    )
    trace = h * float(np.trace(j))
    y = np.zeros(size)
    y[2*(n//2)+1] = 1.0
    exact = y.copy()
    estimate = np.zeros(size)
    oracle = np.zeros(size)
    max_local = 0.0
    elapsed_ms = 0.0
    steps = round(duration/h)
    if not np.isclose(steps*h, duration):
        raise ValueError('duration must be a multiple of dt')
    for _ in range(steps):
        ynext = g @ y
        r0 = (ynext-y)/h - j@y
        r1 = (ynext-y)/h - j@ynext
        # The augmented state evolves for h; its first block yields
        # exp(hJ)*estimate + integral exp((h-s)J)*r(s) ds.
        start = np.concatenate((estimate, r0, (r1-r0)/h))
        t0 = time.perf_counter()
        estimate = expm_multiply(h*action, start, traceA=trace)[:size]
        elapsed_ms += (time.perf_counter()-t0)*1000
        defect = (g-e) @ y
        oracle = e@oracle + defect
        # Independent local defect action (not included in timed run).
        local = expm_multiply(h*augmented,
                             np.concatenate((np.zeros(size), r0, (r1-r0)/h)))[:size]
        local_err = np.linalg.norm(local-defect)/max(np.linalg.norm(defect), 1e-14)
        max_local = max(max_local, local_err)
        y = ynext
        exact = e @ exact
    actual = y-exact
    denom = max(np.linalg.norm(actual), 1e-14)
    return dict(n=n, production_theta=theta, dt_s=h, steps=steps,
                method='scipy-expm-multiply-augmented',
                effectivity=np.linalg.norm(estimate)/denom,
                vector_relative_error=np.linalg.norm(estimate-actual)/denom,
                oracle_relative_error=np.linalg.norm(oracle-actual)/denom,
                elapsed_ms=elapsed_ms, jacobian_products=counts['columns'],
                relative_local_defect_error=max_local)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--duration', type=float, default=600.0)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'results-m7c'/'expm_action_summary.csv')
    args=parser.parse_args()
    rows=[]
    for n,j in sorted(load_operators().items()):
        for theta in (0.65, 1.0):
            for h in (7.5,15.0,30.0,60.0):
                row=run_case(j,n,theta,h,args.duration)
                rows.append(row)
                print(f'{n:2d} {theta:.2f} {h:5.1f} '
                      f'effectivity={row["effectivity"]:.9g} '
                      f'vector_error={row["vector_relative_error"]:.4e} '
                      f'local_error={row["relative_local_defect_error"]:.4e} '
                      f'products={row["jacobian_products"]} '
                      f'time_ms={row["elapsed_ms"]:.1f}', flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(rows)
    print('Max oracle error:',max(float(r['oracle_relative_error']) for r in rows))
    print('Wrote',args.output)

if __name__=='__main__':
    main()
