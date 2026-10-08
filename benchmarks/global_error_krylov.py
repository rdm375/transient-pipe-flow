#!/usr/bin/env python3
"""M7c-6: augmented Arnoldi exponential action for residual-driven error.

Only the validation reference uses a dense expm. Candidate estimators use
J-vector products and exponentials of small Arnoldi Hessenberg matrices.
"""
import argparse
import csv
from time import perf_counter
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import ROOT, load_operators

FIELDS = ('n', 'production_theta', 'dt_s', 'krylov_dimension',
          'effectivity', 'vector_relative_error', 'signed_alignment',
          'j_matvecs', 'arnoldi_seconds', 'projected_expm_seconds',
          'estimator_seconds', 'oracle_relative_error')


def augmented_action(j, error, r0, dr, h, dimension):
    """Integrate e'=J e+r0+(s/h)dr, using one Arnoldi exponential action.

    Augment with a'=0, b'=a, initial (a,b)=(1,0). Then
    e'=J e+r0*a+(dr/h)*b, with b=s.
    """
    n = len(error)
    initial = np.empty(n + 2)
    initial[:n] = error
    initial[n:] = (1., 0.)
    beta = np.linalg.norm(initial)
    if beta == 0:
        return error.copy(), 0, 0., 0.
    maxdim = min(dimension, n+2)
    basis = np.zeros((n+2, maxdim+1))
    hes = np.zeros((maxdim+1, maxdim))
    basis[:, 0] = initial / beta
    count = 0
    start = perf_counter()
    for k in range(maxdim):
        v = basis[:, k]
        w = np.empty_like(v)
        w[:n] = j @ v[:n] + r0*v[n] + (dr/h)*v[n+1]
        w[n] = 0.
        w[n+1] = v[n]
        count += 1
        # Modified Gram-Schmidt with a second pass for nonnormal operators.
        for _ in range(2):
            for i in range(k+1):
                projection = np.dot(basis[:, i], w)
                hes[i, k] += projection
                w -= projection*basis[:, i]
        hes[k+1, k] = np.linalg.norm(w)
        if hes[k+1, k] <= 1e-13 * max(1., np.linalg.norm(hes[:k+1,:k+1])):
            count = k+1
            break
        if k+1 < maxdim:
            basis[:, k+1] = w/hes[k+1,k]
    arnoldi_time = perf_counter()-start
    start = perf_counter()
    small = expm(h*hes[:count,:count])[:, 0]
    answer = beta*(basis[:n,:count] @ small)
    expm_time = perf_counter()-start
    return answer, count, arnoldi_time, expm_time


def run_case(j, n, theta, h, duration, dimensions):
    steps = round(duration/h)
    if abs(steps*h-duration)>1e-9:
        raise ValueError('duration must be divisible by timestep')
    eye = np.eye(len(j))
    lu = lu_factor(eye-theta*h*j)
    g = lu_solve(lu, eye+(1-theta)*h*j)
    exact_step = expm(h*j)  # ORACLE ONLY
    initial = np.zeros(len(j))
    initial[2*(n//2)+1] = 1.
    y = initial.copy()
    reference = initial.copy()
    oracle = np.zeros_like(initial)
    configs = {m: dict(error=np.zeros_like(initial), matvecs=0,
                       arnoldi=0., projected=0.) for m in dimensions}
    for _ in range(steps):
        next_y = g@y
        r0 = (next_y-y)/h-j@y
        dr = -j@(next_y-y)
        for m, c in configs.items():
            c['error'], products, ta, te = augmented_action(
                j, c['error'], r0, dr, h, m)
            c['matvecs'] += products
            c['arnoldi'] += ta
            c['projected'] += te
        oracle = exact_step@oracle + (g-exact_step)@y
        y = next_y
        reference = exact_step@reference
    actual = y-reference
    anorm = np.linalg.norm(actual)
    oracle_relative = np.linalg.norm(oracle-actual)/max(anorm, 1e-14)
    rows = []
    for m, c in configs.items():
        estimate = c['error']
        enorm = np.linalg.norm(estimate)
        rows.append(dict(n=n, production_theta=theta, dt_s=h,
            krylov_dimension=m, effectivity=enorm/anorm,
            vector_relative_error=np.linalg.norm(estimate-actual)/anorm,
            signed_alignment=np.dot(estimate,actual)/(enorm*anorm)
                if enorm*anorm>1e-28 else np.nan,
            j_matvecs=c['matvecs'], arnoldi_seconds=c['arnoldi'],
            projected_expm_seconds=c['projected'],
            estimator_seconds=c['arnoldi']+c['projected'],
            oracle_relative_error=oracle_relative))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dimensions', nargs='+', type=int,
                        default=[8, 16, 24, 32, 48])
    parser.add_argument('--cells', nargs='+', type=int, default=[20,40])
    parser.add_argument('--dt', nargs='+', type=float, default=[30.,60.])
    parser.add_argument('--duration', type=float, default=600.)
    args = parser.parse_args()
    all_rows=[]
    operators=load_operators()
    print('n theta dt krylov effectivity vector-error J-products estimator-ms')
    for n in args.cells:
        j=operators[n]
        for theta in (0.65, 1.):
            for h in args.dt:
                for row in run_case(j,n,theta,h,args.duration,args.dimensions):
                    all_rows.append(row)
                    print(f'{n:2d} {theta:4.2f} {h:5.1f} '
                          f'{row["krylov_dimension"]:6d} '
                          f'{row["effectivity"]:11.6g} '
                          f'{row["vector_relative_error"]:12.6g} '
                          f'{row["j_matvecs"]:10d} '
                          f'{1000*row["estimator_seconds"]:12.3f}')
    max_oracle=max(row['oracle_relative_error'] for row in all_rows)
    print(f'Max oracle relative error: {max_oracle:.3e}')
    if max_oracle>1e-9:
        raise SystemExit('FAIL: oracle transport identity')
    out=ROOT/'results-m7c'/'krylov_summary.csv'
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f'Wrote {out}')
    print('PASS: oracle transport identity')


if __name__=='__main__':
    main()
