#!/usr/bin/env python3
"""M7c-11: offline propagation-weighted Krylov error-budget experiment.

Uses the exact augmented input at every production step and an independent
matrix exponential to measure each candidate action error. This is an oracle
allocation study, NOT a deployable adaptive controller. No production step doubling.

The candidate approximation is evaluated on the exact augmented input, so the
sum of propagated first-block candidate defects is the final error for this
frozen-input experiment. A live adaptive run must additionally account for
feedback of previous approximation errors into later Arnoldi inputs.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from scipy.sparse import bmat, csc_matrix, eye
from scipy.sparse.linalg import expm_multiply

from global_error_transport import ROOT, load_operators
from global_error_adaptive_krylov import initial_state


def candidates(a, v, h, dimensions, size):
    """One incremental Arnoldi basis, all requested dimensions; count m Jv."""
    beta = np.linalg.norm(v)
    if beta == 0:
        return {m: np.zeros(size) for m in dimensions}
    cap = min(max(dimensions), len(v))
    q = np.zeros((len(v), cap + 1))
    hh = np.zeros((cap + 1, cap))
    q[:, 0] = v / beta
    result = {}
    for k in range(cap):
        w = a @ q[:, k]
        for _ in range(2):
            for i in range(k + 1):
                coeff = np.dot(q[:, i], w)
                hh[i, k] += coeff
                w -= coeff * q[:, i]
        m = k + 1
        hnext = np.linalg.norm(w)
        hh[m, k] = hnext
        breakdown = hnext <= 1e-13 * max(1., np.linalg.norm(hh[:m, :m], ord=np.inf))
        if m in dimensions or breakdown:
            result[m] = beta * (q[:, :m] @ expm(h * hh[:m, :m])[:, 0])[:size]
        if breakdown:
            for dim in dimensions:
                if dim >= m:
                    result[dim] = result[m].copy()
            break
        if m < cap:
            q[:, m] = w / hnext
    return {m: result[m] for m in dimensions}


def evaluate_allocation(chosen, defects, propagators, true_norm):
    accumulated = np.zeros_like(defects[0][chosen[0]])
    sum_norm = 0.
    for k, m in enumerate(chosen):
        propagated = propagators[k] @ defects[k][m]
        accumulated += propagated
        sum_norm += np.linalg.norm(propagated)
    return np.linalg.norm(accumulated) / true_norm, sum_norm / true_norm


def allocate(defects, propagators, dimensions, true_norm, target, method):
    steps = len(defects)
    m_min, m_max = dimensions[0], dimensions[-1]
    if method == 'uniform':
        # Equal dimension everywhere: smallest choice satisfying global target.
        for m in dimensions:
            chosen = [m] * steps
            error, bound = evaluate_allocation(chosen, defects, propagators, true_norm)
            if bound <= target:
                return chosen, error, bound
        return [m_max] * steps, *evaluate_allocation([m_max] * steps, defects, propagators, true_norm)

    # Greedy cost-aware oracle: reduce the triangle-inequality bound per Jv.
    # 'local' uses ||delta_k||; 'weighted' uses ||E^(N-k-1) delta_k||.
    chosen = [m_min] * steps
    def score(k, m):
        vector = defects[k][m]
        if method == 'weighted':
            vector = propagators[k] @ vector
        return np.linalg.norm(vector)
    current = [score(k, m_min) for k in range(steps)]
    while sum(current) > target * true_norm:
        best = None
        for k in range(steps):
            idx = dimensions.index(chosen[k])
            if idx + 1 >= len(dimensions):
                continue
            nxt = dimensions[idx + 1]
            next_score = score(k, nxt)
            saving = current[k] - next_score
            efficiency = saving / (nxt - chosen[k])
            if best is None or efficiency > best[0]:
                best = (efficiency, k, nxt, next_score)
        if best is None or best[0] <= 0:
            break
        _, k, nxt, new_score = best
        chosen[k] = nxt
        current[k] = new_score
    error, bound = evaluate_allocation(chosen, defects, propagators, true_norm)
    return chosen, error, bound


def run_case(j, n, theta, dt, duration, probe, seed, dimensions, targets):
    size = len(j)
    ident = np.eye(size)
    g = lu_solve(lu_factor(ident - theta * dt * j), ident + (1 - theta) * dt * j)
    e = expm(dt * j)
    zero = csc_matrix((size, size))
    a = bmat([[csc_matrix(j), eye(size, format='csc'), zero],
              [zero, zero, eye(size, format='csc')],
              [zero, zero, zero]], format='csc')
    steps = round(duration / dt)
    if not np.isclose(steps * dt, duration):
        raise ValueError('duration must be divisible by dt')
    y = initial_state(j, n, probe, seed)
    exact = y.copy()
    estimator = np.zeros(size)
    defects = []
    propagators = []
    for k in range(steps):
        ynext = g @ y
        r0 = (ynext - y) / dt - j @ y
        r1 = (ynext - y) / dt - j @ ynext
        v = np.concatenate((estimator, r0, (r1 - r0) / dt))
        oracle = expm_multiply(dt * a, v)[:size]
        approx = candidates(a, v, dt, dimensions, size)
        defects.append({m: approx[m] - oracle for m in dimensions})
        estimator = oracle
        exact = e @ exact
        y = ynext
    actual = y - exact
    true_norm = np.linalg.norm(actual)
    oracle_discrepancy = np.linalg.norm(estimator - actual)
    if true_norm <= 1e-14:
        raise ValueError('true global error too small for meaningful relative budgeting')
    for k in range(steps):
        propagators.append(np.linalg.matrix_power(e, steps - 1 - k))
    rows = []
    for target in targets:
        for method in ('uniform', 'local', 'weighted'):
            selected, actual_relative, bound_relative = allocate(
                defects, propagators, dimensions, true_norm, target, method)
            row = dict(n=n, theta=theta, dt_s=dt, probe=probe,
                       target_relative=target, policy=method,
                       frozen_global_relative_error=actual_relative,
                       triangle_bound_relative=bound_relative,
                       cost_jv=sum(selected), mean_dimension=np.mean(selected),
                       max_dimension=max(selected), dimensions=';'.join(map(str, selected)),
                       true_global_norm=true_norm,
                       exact_recurrence_error=oracle_discrepancy,
                       target_met_bound=int(bound_relative <= target),
                       target_met_actual=int(actual_relative <= target))
            rows.append(row)
            print(f'{n} theta={theta:g} dt={dt:g} {probe:7s} target={target:.0e} '
                  f'{method:8s} Jv={row["cost_jv"]:4d} '
                  f'actual={actual_relative:.4g} bound={bound_relative:.4g} '
                  f'bound_met={row["target_met_bound"]} '
                  f'dims={row["dimensions"]}', flush=True)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, nargs='+', default=[40])
    p.add_argument('--theta', type=float, nargs='+', default=[.65])
    p.add_argument('--dt', type=float, nargs='+', default=[60.])
    p.add_argument('--probes', nargs='+', default=['fast', 'slow', 'impulse'],
                   choices=['fast', 'slow', 'impulse', 'random'])
    p.add_argument('--targets', type=float, nargs='+', default=[1e-2, 1e-4])
    p.add_argument('--max-dim', type=int, default=32)
    p.add_argument('--checkpoint', type=int, default=4)
    p.add_argument('--duration', type=float, default=600.)
    p.add_argument('--seed', type=int, default=1729)
    p.add_argument('--output', type=Path, default=ROOT / 'results-m7c' / 'global_budgeting.csv')
    args = p.parse_args()
    if args.max_dim < 8 or args.checkpoint < 1 or any(x <= 0 for x in args.targets):
        p.error('require max-dim>=8, checkpoint>=1, and positive targets')
    dims = list(range(max(4, args.checkpoint), args.max_dim + 1, args.checkpoint))
    if dims[-1] != args.max_dim:
        dims.append(args.max_dim)
    operators = load_operators()
    rows = []
    for n in args.n:
        for theta in args.theta:
            for dt in args.dt:
                for probe in args.probes:
                    rows.extend(run_case(operators[n], n, theta, dt, args.duration,
                                         probe, args.seed, dims, args.targets))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print('Wrote', args.output)
    print('ORACLE/FROZEN-INPUT STUDY: Jv is sum of selected basis dimensions,')
    print('excluding oracle evaluations and candidate exploration.')
    print('If bound_met=0, the dimension cap prevented meeting the target.')


if __name__ == '__main__':
    main()
