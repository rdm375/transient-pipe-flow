#!/usr/bin/env python3
"""M7c-8: incremental Arnoldi, checkpoint diagnostics and independent oracle.

Experimental controller, NOT a certified error bound. No step doubling of the
production solver. The dense exponential is used only for oracle validation.
"""
import argparse
import csv
import time
from pathlib import Path
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import ROOT, load_operators

FIELDS = ('n','production_theta','dt_s','probe','requested_rtol','effectivity',
          'vector_relative_error','jacobian_products','accepted_actions',
          'rejected_actions','max_krylov_dimension','mean_accepted_dimension',
          'elapsed_ms','oracle_relative_error','max_indicator_ratio')


def incremental_action(matvec, v, t, rtol, dimensions, counter, stats,
                       depth=0, max_depth=16):
    """Single expanding Arnoldi basis; residual indicator is NOT a bound."""
    beta = np.linalg.norm(v)
    if beta == 0:
        stats['accepted'] += 1
        return v.copy()
    n = len(v)
    cap = min(dimensions[-1], n)
    q = np.zeros((n, cap+1))
    hh = np.zeros((cap+1, cap))
    q[:,0] = v / beta
    previous = None
    last = None
    for k in range(cap):
        w = matvec(q[:,k]); counter[0] += 1
        for _ in range(2):
            for i in range(k+1):
                c = np.dot(q[:,i],w)
                hh[i,k] += c
                w -= c*q[:,i]
        hnext = np.linalg.norm(w)
        hh[k+1,k] = hnext
        m = k+1
        breakdown = hnext <= 1e-13 * max(1.0, np.linalg.norm(hh[:m,:m],ord=np.inf))
        if not breakdown and m < cap:
            q[:,m] = w/hnext
        if m not in dimensions and not breakdown and m != cap:
            continue
        hm = hh[:m,:m]
        # One (m+1)-dimensional exponential provides the approximation and
        # signed integrated Arnoldi residual coefficient (heuristic only).
        aug = np.zeros((m+1,m+1))
        aug[:m,:m] = hm
        aug[m,m-1] = 0.0 if breakdown else hnext
        col = beta*expm(t*aug)[:,0]
        current = q[:,:m] @ col[:m]
        indicator = abs(col[m])
        scale = max(np.linalg.norm(current), np.linalg.norm(v), 1e-30)
        ratio = indicator/scale
        stats['max_indicator_ratio'] = max(stats['max_indicator_ratio'], ratio)
        # Difference with previous checkpoint is an additional diagnostic;
        # cannot be interpreted as a rigorous convergence certificate.
        difference = (np.linalg.norm(current-previous)/scale
                      if previous is not None else np.inf)
        stats['max_dimension'] = max(stats['max_dimension'],m)
        last = current
        if breakdown or (previous is not None and
                         ratio <= rtol and difference <= rtol):
            stats['accepted'] += 1
            stats['accepted_dimensions'].append(m)
            return current
        previous = current
    stats['rejected'] += 1
    if depth >= max_depth:
        raise RuntimeError('Krylov action failed to converge: increase --max-dim or --max-depth')
    # This subdivides ONLY the error-estimator action, not the production step.
    mid = incremental_action(matvec,v,t/2,rtol/2,dimensions,counter,stats,
                             depth+1,max_depth)
    return incremental_action(matvec,mid,t/2,rtol/2,dimensions,counter,stats,
                              depth+1,max_depth)


def initial_state(j,n,probe,seed):
    size=len(j)
    if probe == 'impulse':
        y=np.zeros(size); y[2*(n//2)+1]=1.; return y
    if probe == 'random':
        rng=np.random.default_rng(seed+n)
        y=rng.standard_normal(size); return y/np.linalg.norm(y)
    if probe in ('slow','fast'):
        eigvals, eigvecs=np.linalg.eig(j)
        # Select spectral modes using real decay rate, normalize real part.
        order=np.argsort(np.real(eigvals))
        for idx in (order[::-1] if probe=='slow' else order):
            vec=np.real(eigvecs[:,idx])
            if np.linalg.norm(vec)>1e-12:
                return vec/np.linalg.norm(vec)
    raise ValueError(probe)


def run_case(j,n,theta,dt,duration,tol,dimensions,probe,seed,max_depth):
    size=len(j); ident=np.eye(size)
    g=lu_solve(lu_factor(ident-theta*dt*j),ident+(1-theta)*dt*j)
    e=expm(dt*j)  # oracle only
    steps=round(duration/dt)
    if not np.isclose(steps*dt,duration):
        raise ValueError('duration must be a multiple of dt')
    y=initial_state(j,n,probe,seed)
    exact=y.copy(); oracle=np.zeros(size); estimate=np.zeros(size)
    count=[0]
    stats=dict(accepted=0,rejected=0,max_dimension=0,
               accepted_dimensions=[],max_indicator_ratio=0.)
    def product(x):
        z,v,w=np.split(x,3)
        return np.concatenate((j@z+v,w,np.zeros_like(w)))
    t0=time.perf_counter()
    for _ in range(steps):
        ynext=g@y
        r0=(ynext-y)/dt-j@y
        r1=(ynext-y)/dt-j@ynext
        start=np.concatenate((estimate,r0,(r1-r0)/dt))
        estimate=incremental_action(product,start,dt,tol,dimensions,count,
                                    stats,max_depth=max_depth)[:size]
        oracle=e@oracle+(g-e)@y
        exact=e@exact
        y=ynext
    ms=1000*(time.perf_counter()-t0)
    actual=y-exact
    denom=max(np.linalg.norm(actual),1e-14)
    return dict(n=n,production_theta=theta,dt_s=dt,probe=probe,
                requested_rtol=tol,effectivity=np.linalg.norm(estimate)/denom,
                vector_relative_error=np.linalg.norm(estimate-actual)/denom,
                jacobian_products=count[0],accepted_actions=stats['accepted'],
                rejected_actions=stats['rejected'],
                max_krylov_dimension=stats['max_dimension'],
                mean_accepted_dimension=np.mean(stats['accepted_dimensions']),
                elapsed_ms=ms,oracle_relative_error=np.linalg.norm(oracle-actual)/denom,
                max_indicator_ratio=stats['max_indicator_ratio'])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--n',type=int,nargs='+',default=[20,40])
    p.add_argument('--dt',type=float,nargs='+',default=[30.,60.])
    p.add_argument('--theta',type=float,nargs='+',default=[.65,1.])
    p.add_argument('--tol',type=float,nargs='+',default=[1e-2,1e-4,1e-6])
    p.add_argument('--probes',nargs='+',default=['impulse','random','slow','fast'],
                   choices=['impulse','random','slow','fast'])
    p.add_argument('--seed',type=int,default=1729)
    p.add_argument('--max-dim',type=int,default=32)
    p.add_argument('--checkpoint',type=int,default=4)
    p.add_argument('--max-depth',type=int,default=16)
    p.add_argument('--duration',type=float,default=600.)
    p.add_argument('--output',type=Path,
                   default=ROOT/'results-m7c'/'adaptive_krylov_summary.csv')
    args=p.parse_args()
    if args.max_dim<8 or args.checkpoint<1:
        p.error('--max-dim must be >=8 and --checkpoint must be >=1')
    dims=list(range(max(4,args.checkpoint),args.max_dim+1,args.checkpoint))
    if args.max_dim not in dims: dims.append(args.max_dim)
    operators=load_operators(); rows=[]
    for n in args.n:
        for theta in args.theta:
            for dt in args.dt:
                for probe in args.probes:
                    for tol in args.tol:
                        row=run_case(operators[n],n,theta,dt,args.duration,
                                     tol,dims,probe,args.seed,args.max_depth)
                        rows.append(row)
                        print(f'{n:2d} {theta:.2f} {dt:5g} {probe:7s} '
                              f'tol={tol:.0e} err={row["vector_relative_error"]:.4g} '
                              f'Jv={row["jacobian_products"]} '
                              f'mean_m={row["mean_accepted_dimension"]:.1f} '
                              f'reject={row["rejected_actions"]} '
                              f'ms={row["elapsed_ms"]:.1f}',flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader(); writer.writerows(rows)
    print('Wrote',args.output)
    print('Max oracle error:',max(r['oracle_relative_error'] for r in rows))

if __name__=='__main__':
    main()
