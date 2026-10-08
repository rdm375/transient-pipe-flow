#!/usr/bin/env python3
"""M7c-10: compare Krylov stopping rules on identical operator/probe cases.

Controllers:
  baseline: M7c-8 full-state residual + checkpoint test.
  output: first-block checkpoint test (heuristic, no oracle in decisions).
  propagated_oracle: true propagated first-block Arnoldi error, computed
      with an independent exponential action. Diagnostic LOWER BOUND ON COST,
      not a deployable algorithm or certified a posteriori bound.

The oracle criterion is equivalent to evaluating the norm of the full
propagated Arnoldi residual integral, not an unpropagated residual scalar.
No step doubling of the production integrator.
"""
import argparse
import csv
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from scipy.sparse import bmat, csc_matrix, eye
from scipy.sparse.linalg import expm_multiply

from global_error_transport import ROOT, load_operators
from global_error_adaptive_krylov import initial_state

CONTROLLERS = ('baseline', 'output', 'propagated_oracle')
FIELDS = ('n','production_theta','dt_s','probe','requested_rtol','controller',
          'true_global_norm','absolute_global_error','relative_global_error',
          'effectivity','jacobian_products','accepted_actions','rejected_actions',
          'mean_dimension','max_dimension','worst_accepted_first_action_error',
          'worst_accepted_first_action_relative_error','false_acceptances',
          'oracle_action_calls','elapsed_ms')


def arnoldi_controlled(a, v, t, tol, dims, controller, counts, size,
                       depth=0, max_depth=16, oracle=None):
    beta = np.linalg.norm(v)
    if beta == 0:
        counts['accepted'] += 1
        counts['dimensions'].append(0)
        return v.copy()
    cap = min(dims[-1], len(v))
    q = np.zeros((len(v), cap+1))
    hh = np.zeros((cap+1, cap))
    q[:, 0] = v/beta
    previous = None
    # Reference only for diagnostics. For the oracle controller it is also
    # the acceptance criterion; never used to form the returned approximation.
    if oracle is None:
        oracle = expm_multiply(t*a, v)
        counts['oracle_calls'] += 1
    ref_first = oracle[:size]
    for k in range(cap):
        w = a @ q[:, k]
        counts['jv'] += 1
        for _ in range(2):
            for i in range(k+1):
                c = np.dot(q[:, i], w)
                hh[i, k] += c
                w -= c*q[:, i]
        m = k+1
        hnext = np.linalg.norm(w)
        hh[m, k] = hnext
        breakdown = hnext <= 1e-13 * max(1., np.linalg.norm(hh[:m,:m],ord=np.inf))
        if not breakdown and m < cap:
            q[:, m] = w/hnext
        if m not in dims and not breakdown and m != cap:
            continue
        small = np.zeros((m+1,m+1))
        small[:m,:m] = hh[:m,:m]
        small[m,m-1] = 0. if breakdown else hnext
        col = beta*expm(t*small)[:,0]
        current = q[:,:m] @ col[:m]
        indicator = abs(col[m])
        first_error = np.linalg.norm(current[:size]-ref_first)
        first_scale = max(np.linalg.norm(current[:size]),
                          np.linalg.norm(v[:size]), 1e-30)
        oracle_scale = max(np.linalg.norm(ref_first),
                           np.linalg.norm(v[:size]), 1e-30)
        if controller == 'baseline':
            scale = max(np.linalg.norm(current), beta, 1e-30)
            delta = np.inf if previous is None else np.linalg.norm(current-previous)
            accept = indicator <= tol*scale and delta <= tol*scale
        elif controller == 'output':
            delta = np.inf if previous is None else np.linalg.norm(current[:size]-previous[:size])
            # The residual scalar is NOT an output bound. Require only
            # first-block checkpoint agreement, which remains heuristic.
            accept = delta <= tol*first_scale
        else:
            # Exact propagated Arnoldi error (diagnostic oracle only).
            accept = first_error <= tol*oracle_scale
        accept = breakdown or accept
        if accept:
            counts['accepted'] += 1
            counts['dimensions'].append(m)
            counts['worst_abs'] = max(counts['worst_abs'], first_error)
            counts['worst_rel'] = max(counts['worst_rel'], first_error/oracle_scale)
            if first_error > tol*oracle_scale:
                counts['false_accept'] += 1
            return current
        previous = current
    counts['rejected'] += 1
    if depth >= max_depth:
        raise RuntimeError('Krylov controller failed: raise --max-dim or --max-depth')
    mid = arnoldi_controlled(a,v,t/2,tol/2,dims,controller,counts,size,
                             depth+1,max_depth)
    return arnoldi_controlled(a,mid,t/2,tol/2,dims,controller,counts,size,
                              depth+1,max_depth)


def run_case(j,n,theta,dt,duration,tol,dims,probe,seed,controller,max_depth):
    size = len(j)
    ident = np.eye(size)
    g = lu_solve(lu_factor(ident-theta*dt*j), ident+(1-theta)*dt*j)
    e = expm(dt*j)  # Independent reference only.
    z = csc_matrix((size,size))
    a = bmat([[csc_matrix(j),eye(size,format='csc'),z],
              [z,z,eye(size,format='csc')],[z,z,z]],format='csc')
    steps = round(duration/dt)
    if not np.isclose(steps*dt,duration):
        raise ValueError('duration must be a multiple of dt')
    y = initial_state(j,n,probe,seed)
    exact = y.copy()
    estimate = np.zeros(size)
    counts = dict(jv=0,accepted=0,rejected=0,dimensions=[],worst_abs=0.,
                  worst_rel=0.,false_accept=0,oracle_calls=0)
    t0 = time.perf_counter()
    for _ in range(steps):
        ynext = g@y
        r0 = (ynext-y)/dt-j@y
        r1 = (ynext-y)/dt-j@ynext
        start = np.concatenate((estimate,r0,(r1-r0)/dt))
        estimate = arnoldi_controlled(a,start,dt,tol,dims,controller,counts,size,
                                      max_depth=max_depth)[:size]
        exact = e@exact
        y = ynext
    elapsed = 1000*(time.perf_counter()-t0)
    actual = y-exact
    true_norm = np.linalg.norm(actual)
    denom = max(true_norm,1e-14)
    return dict(n=n,production_theta=theta,dt_s=dt,probe=probe,
                requested_rtol=tol,controller=controller,true_global_norm=true_norm,
                absolute_global_error=np.linalg.norm(estimate-actual),
                relative_global_error=np.linalg.norm(estimate-actual)/denom,
                effectivity=np.linalg.norm(estimate)/denom,
                jacobian_products=counts['jv'],accepted_actions=counts['accepted'],
                rejected_actions=counts['rejected'],
                mean_dimension=np.mean(counts['dimensions']) if counts['dimensions'] else 0,
                max_dimension=max(counts['dimensions'],default=0),
                worst_accepted_first_action_error=counts['worst_abs'],
                worst_accepted_first_action_relative_error=counts['worst_rel'],
                false_acceptances=counts['false_accept'],
                oracle_action_calls=counts['oracle_calls'],elapsed_ms=elapsed)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,nargs='+',default=[40])
    p.add_argument('--theta',type=float,nargs='+',default=[.65])
    p.add_argument('--dt',type=float,nargs='+',default=[60.])
    p.add_argument('--probes',nargs='+',default=['fast','slow','impulse'],
                   choices=['impulse','random','slow','fast'])
    p.add_argument('--tol',type=float,nargs='+',default=[1e-2,1e-4])
    p.add_argument('--controllers',nargs='+',choices=CONTROLLERS,default=list(CONTROLLERS))
    p.add_argument('--seed',type=int,default=1729)
    p.add_argument('--max-dim',type=int,default=32)
    p.add_argument('--checkpoint',type=int,default=4)
    p.add_argument('--max-depth',type=int,default=16)
    p.add_argument('--duration',type=float,default=600.)
    p.add_argument('--output',type=Path,default=ROOT/'results-m7c'/'controller_comparison.csv')
    args=p.parse_args()
    if args.max_dim<8 or args.checkpoint<1:
        p.error('max-dim must be >=8 and checkpoint >=1')
    dims=list(range(max(4,args.checkpoint),args.max_dim+1,args.checkpoint))
    if args.max_dim not in dims: dims.append(args.max_dim)
    operators=load_operators()
    rows=[]
    for n in args.n:
        for theta in args.theta:
            for dt in args.dt:
                for probe in args.probes:
                    for tol in args.tol:
                        for controller in args.controllers:
                            row=run_case(operators[n],n,theta,dt,args.duration,tol,dims,
                                         probe,args.seed,controller,args.max_depth)
                            rows.append(row)
                            print(f'{n} theta={theta:g} dt={dt:g} {probe:7s} '
                                  f'tol={tol:.0e} {controller:18s} '
                                  f'err={row["relative_global_error"]:.5g} '
                                  f'Jv={row["jacobian_products"]:4d} '
                                  f'mean_m={row["mean_dimension"]:.1f} '
                                  f'false={row["false_acceptances"]} '
                                  f'reject={row["rejected_actions"]} '
                                  f'ms={row["elapsed_ms"]:.1f}',flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(rows)
    print('Wrote',args.output)
    print('Note: elapsed_ms INCLUDES expensive oracle expm_multiply calls for all controllers.')
    print('Only Jacobian products in the incremental Arnoldi algorithm are counted in Jv.')

if __name__=='__main__':
    main()
