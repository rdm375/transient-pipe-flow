#!/usr/bin/env python3
"""M7c-12: exact discrete Pareto audit of frozen-input Krylov allocations.

For each production step, construct one Arnoldi basis and compare every
checkpoint with an independent exponential-action oracle. Find the EXACT
minimum sum of checkpoint dimensions satisfying the triangle bound by
multiple-choice dynamic programming indexed by integer Jv cost.

The experiment is offline/oracle-assisted. Costs exclude discovering the
candidate table, evaluating oracles, and feedback in a live estimator.
No production step doubling.
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
from global_error_budgeting import candidates, evaluate_allocation, allocate


def pareto_dp(costs, errors):
    """Exact minimal sum of nonnegative errors at each achievable Jv cost."""
    stages = [{0: (0., ())}]
    for row in errors:
        prev = stages[-1]
        curr = {}
        for cost0, (error0, path0) in prev.items():
            for m, error in zip(costs, row):
                cost = cost0 + m
                value = error0 + float(error)
                if cost not in curr or value < curr[cost][0]:
                    curr[cost] = (value, path0 + (m,))
        stages.append(curr)
    return stages[-1]


def run_case(j, n, theta, dt, duration, probe, seed, dimensions, targets):
    size = len(j)
    ident = np.eye(size)
    g = lu_solve(lu_factor(ident - theta * dt * j), ident + (1 - theta) * dt * j)
    e = expm(dt * j)
    z = csc_matrix((size, size))
    a = bmat([[csc_matrix(j), eye(size, format='csc'), z],
              [z, z, eye(size, format='csc')], [z, z, z]], format='csc')
    steps = round(duration / dt)
    if not np.isclose(steps * dt, duration):
        raise ValueError('duration must be divisible by dt')
    y = initial_state(j, n, probe, seed)
    exact = y.copy()
    estimator = np.zeros(size)
    defects = []
    local_norms = []
    propagators = []
    for k in range(steps):
        ynext = g @ y
        r0 = (ynext - y) / dt - j @ y
        r1 = (ynext - y) / dt - j @ ynext
        v = np.concatenate((estimator, r0, (r1 - r0) / dt))
        oracle = expm_multiply(dt * a, v)[:size]
        approx = candidates(a, v, dt, dimensions, size)
        defects.append({m: approx[m] - oracle for m in dimensions})
        local_norms.append([np.linalg.norm(defects[-1][m]) for m in dimensions])
        estimator = oracle
        exact = e @ exact
        y = ynext
    actual = y - exact
    true_norm = np.linalg.norm(actual)
    oracle_gap = np.linalg.norm(estimator - actual)
    if true_norm <= 1e-14:
        raise ValueError('true global error too small for relative budgeting')
    for k in range(steps):
        propagators.append(np.linalg.matrix_power(e, steps - 1 - k))
    weighted = [[np.linalg.norm(propagators[k] @ defects[k][m])
                 for m in dimensions] for k in range(steps)]
    dp = pareto_dp(dimensions, weighted)
    frontier = []
    best = np.inf
    for cost in sorted(dp):
        bound, chosen = dp[cost]
        if bound < best - 1e-15 * max(1., best if np.isfinite(best) else 1.):
            actual_rel, bound_rel = evaluate_allocation(chosen, defects, propagators, true_norm)
            frontier.append(dict(n=n, theta=theta, dt_s=dt, probe=probe,
                                 cost_jv=cost, triangle_bound_relative=bound_rel,
                                 frozen_actual_relative=actual_rel,
                                 dimensions=';'.join(map(str, chosen))))
            best = bound
    rows = []
    for target in targets:
        feasible = [(cost, val) for cost, val in dp.items()
                    if val[0] <= target * true_norm]
        if feasible:
            cost_opt, (_, chosen_opt) = min(feasible, key=lambda item: item[0])
        else:
            cost_opt, (_, chosen_opt) = min(dp.items(), key=lambda item: item[1][0])
        for policy in ('uniform', 'local_greedy', 'weighted_greedy', 'optimal_dp'):
            if policy == 'optimal_dp':
                chosen = chosen_opt
            else:
                kind = {'uniform': 'uniform', 'local_greedy': 'local',
                        'weighted_greedy': 'weighted'}[policy]
                chosen, _, _ = allocate(defects, propagators, dimensions,
                                        true_norm, target, kind)
            actual_rel, bound_rel = evaluate_allocation(chosen, defects, propagators, true_norm)
            row = dict(n=n, theta=theta, dt_s=dt, probe=probe, target_relative=target,
                       policy=policy, cost_jv=sum(chosen), frozen_actual_relative=actual_rel,
                       triangle_bound_relative=bound_rel,
                       bound_met=int(bound_rel <= target * (1 + 1e-10)),
                       actual_met=int(actual_rel <= target),
                       dimensions=';'.join(map(str, chosen)),
                       true_global_norm=true_norm, oracle_recurrence_gap=oracle_gap,
                       optimal_feasible=int(bool(feasible)))
            rows.append(row)
            print(f'{n} theta={theta:g} dt={dt:g} {probe:7s} target={target:.0e} '
                  f'{policy:16s} Jv={row["cost_jv"]:3d} '
                  f'actual={actual_rel:.5g} bound={bound_rel:.5g} '
                  f'bound_met={row["bound_met"]} dims={row["dimensions"]}', flush=True)
    diagnostics = []
    for k in range(steps):
        for idx, m in enumerate(dimensions):
            diagnostics.append(dict(n=n, theta=theta, dt_s=dt, probe=probe,
                                    step=k + 1, dimension=m,
                                    local_error_norm=local_norms[k][idx],
                                    propagated_error_norm=weighted[k][idx],
                                    local_error_relative=local_norms[k][idx] / true_norm,
                                    propagated_error_relative=weighted[k][idx] / true_norm))
    nonmono_local = sum(local_norms[k][i+1] > local_norms[k][i] * (1 + 1e-10)
                        for k in range(steps) for i in range(len(dimensions)-1))
    nonmono_weighted = sum(weighted[k][i+1] > weighted[k][i] * (1 + 1e-10)
                           for k in range(steps) for i in range(len(dimensions)-1))
    print(f'  audit: {len(dp)} reachable costs; {len(frontier)} Pareto points; '
          f'nonmonotone local={nonmono_local} weighted={nonmono_weighted}; '
          f'oracle recurrence gap={oracle_gap:.3e}', flush=True)
    return rows, diagnostics, frontier


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('Wrote', path)


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
    p.add_argument('--output-dir', type=Path, default=ROOT / 'results-m7c')
    args = p.parse_args()
    if args.max_dim < 4 or args.checkpoint < 1 or any(x <= 0 for x in args.targets):
        p.error('require max-dim>=4, checkpoint>=1, positive targets')
    dims = list(range(4, args.max_dim + 1, args.checkpoint))
    if dims[-1] != args.max_dim:
        dims.append(args.max_dim)
    dims = sorted(set(dims))
    operators = load_operators()
    rows, candidates_rows, frontier_rows = [], [], []
    for n in args.n:
        for theta in args.theta:
            for dt in args.dt:
                for probe in args.probes:
                    r, c, f = run_case(operators[n], n, theta, dt, args.duration,
                                       probe, args.seed, dims, args.targets)
                    rows.extend(r)
                    candidates_rows.extend(c)
                    frontier_rows.extend(f)
    write_csv(args.output_dir / 'allocation_audit.csv', rows)
    write_csv(args.output_dir / 'allocation_candidates.csv', candidates_rows)
    write_csv(args.output_dir / 'allocation_pareto.csv', frontier_rows)
    print('DP is exact for the listed candidate dimensions and additive triangle bound.')
    print('ORACLE/FROZEN-INPUT STUDY: Jv excludes candidate exploration and oracle costs.')


if __name__ == '__main__':
    main()
