#!/usr/bin/env python3
"""M7c-7: tolerance-controlled Arnoldi actions for reconstructed residual transport.

Experimental convergence controller, NOT a certified error bound. The exact
exponential is used only for independent validation of the final error vector.
"""
import argparse
import csv
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import ROOT, load_operators

FIELDS = ('n', 'production_theta', 'dt_s', 'requested_rtol', 'effectivity',
          'vector_relative_error', 'jacobian_products', 'accepted_actions',
          'rejected_actions', 'max_krylov_dimension', 'elapsed_ms',
          'oracle_relative_error')


def arnoldi_action(a, v, t, mmax, counter):
    """Approximate exp(t*A)v; report residual-integral error indicator.

    Arnoldi residual estimate is heuristic for this nonnormal operator.
    """
    beta = np.linalg.norm(v)
    if beta == 0:
        return v.copy(), 0.0, 0
    dim = len(v)
    mmax = min(mmax, dim)
    q = np.zeros((dim, mmax + 1))
    hess = np.zeros((mmax + 1, mmax))
    q[:, 0] = v / beta
    m = mmax
    for k in range(mmax):
        w = a @ q[:, k]
        counter[0] += 1
        # Two-pass modified Gram-Schmidt to reduce orthogonality loss.
        for _ in range(2):
            for i in range(k + 1):
                z = np.dot(q[:, i], w)
                hess[i, k] += z
                w -= z * q[:, i]
        hess[k + 1, k] = np.linalg.norm(w)
        if hess[k + 1, k] < 1e-14:
            m = k + 1
            break
        if k + 1 < mmax:
            q[:, k + 1] = w / hess[k + 1, k]
    hm = hess[:m, :m]
    # Integrate the Arnoldi residual indicator using an augmented matrix:
    # integral_0^t |h_{m+1,m} e_m^T exp(sH)e1| ds is bounded here
    # heuristically by the norm of an integrated signed residual. This is
    # NOT a rigorous bound; validation below is mandatory.
    aug = np.zeros((m + 1, m + 1))
    aug[:m, :m] = hm
    aug[m, m - 1] = hess[m, m - 1] if m < len(hess) else 0.0
    vec = expm(t * aug)[:, 0] * beta
    return q[:, :m] @ vec[:m], abs(vec[m]), m


def controlled_action(a, v, t, rtol, mmax, counter, stats,
                      depth=0, max_depth=18):
    """Compare two Krylov dimensions, split the interval if they disagree.

    Disagreement is a convergence *indicator*, not a guaranteed error bound.
    """
    lo, _, ml = arnoldi_action(a, v, t, max(4, mmax // 2), counter)
    hi, indicator, mh = arnoldi_action(a, v, t, mmax, counter)
    disagreement = np.linalg.norm(hi - lo)
    scale = max(np.linalg.norm(hi), np.linalg.norm(v), 1e-30)
    # Both independent dimension comparison and Arnoldi residual indicator.
    accept = (disagreement <= rtol * scale and
              indicator <= rtol * scale)
    stats['max_dimension'] = max(stats['max_dimension'], mh)
    if accept:
        stats['accepted'] += 1
        return hi
    stats['rejected'] += 1
    if depth >= max_depth:
        raise RuntimeError('Krylov controller did not converge: increase --max-dim')
    middle = controlled_action(a, v, t/2, rtol/2, mmax,
                               counter, stats, depth+1, max_depth)
    return controlled_action(a, middle, t/2, rtol/2, mmax,
                             counter, stats, depth+1, max_depth)


def run_case(j, n, theta, dt, duration, tol, max_dim):
    size = j.shape[0]
    eye = np.eye(size)
    g = lu_solve(lu_factor(eye-theta*dt*j), eye+(1-theta)*dt*j)
    e = expm(dt*j)  # oracle only
    # A maps [z,v,w] to [Jz+v,w,0].
    def augmented_product(x):
        z, v, w = np.split(x, 3)
        return np.concatenate((j @ z + v, w, np.zeros_like(w)))
    from scipy.sparse.linalg import LinearOperator
    a = LinearOperator((3*size, 3*size), matvec=augmented_product,
                       dtype=float)
    steps = round(duration/dt)
    if not np.isclose(steps*dt, duration):
        raise ValueError('duration must be a multiple of dt')
    y = np.zeros(size)
    y[2*(n//2)+1] = 1.0
    exact = y.copy()
    oracle = np.zeros(size)
    estimate = np.zeros(size)
    count = [0]
    stats = dict(accepted=0, rejected=0, max_dimension=0)
    start_time = time.perf_counter()
    for _ in range(steps):
        ynext = g @ y
        r0 = (ynext-y)/dt-j@y
        r1 = (ynext-y)/dt-j@ynext
        start = np.concatenate((estimate, r0, (r1-r0)/dt))
        estimate = controlled_action(a, start, dt, tol, max_dim,
                                     count, stats)[:size]
        oracle = e @ oracle + (g-e) @ y
        exact = e @ exact
        y = ynext
    elapsed_ms = (time.perf_counter()-start_time)*1000
    actual = y-exact
    denom = max(np.linalg.norm(actual), 1e-14)
    return dict(n=n, production_theta=theta, dt_s=dt,
                requested_rtol=tol,
                effectivity=np.linalg.norm(estimate)/denom,
                vector_relative_error=np.linalg.norm(estimate-actual)/denom,
                jacobian_products=count[0],
                accepted_actions=stats['accepted'],
                rejected_actions=stats['rejected'],
                max_krylov_dimension=stats['max_dimension'],
                elapsed_ms=elapsed_ms,
                oracle_relative_error=np.linalg.norm(oracle-actual)/denom)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, nargs='+', default=[20,40])
    parser.add_argument('--dt', type=float, nargs='+', default=[30.,60.])
    parser.add_argument('--theta', type=float, nargs='+', default=[0.65,1.0])
    parser.add_argument('--tol', type=float, nargs='+', default=[1e-2,1e-4,1e-6])
    parser.add_argument('--max-dim', type=int, default=32)
    parser.add_argument('--duration', type=float, default=600.)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'results-m7c'/'tolerance_action_summary.csv')
    args=parser.parse_args()
    operators=load_operators()
    rows=[]
    for n in args.n:
        for theta in args.theta:
            for dt in args.dt:
                for tol in args.tol:
                    row=run_case(operators[n],n,theta,dt,args.duration,tol,args.max_dim)
                    rows.append(row)
                    print(f'{n:2d} {theta:.2f} {dt:5g} tol={tol:.0e} '
                          f'effectivity={row["effectivity"]:.7g} '
                          f'vector_error={row["vector_relative_error"]:.4g} '
                          f'Jv={row["jacobian_products"]} '
                          f'reject={row["rejected_actions"]} '
                          f'ms={row["elapsed_ms"]:.1f}',flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(rows)
    print('Wrote',args.output)
    print('Max oracle error:',max(r['oracle_relative_error'] for r in rows))

if __name__=='__main__':
    main()
