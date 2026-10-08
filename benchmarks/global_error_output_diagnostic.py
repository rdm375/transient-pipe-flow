#!/usr/bin/env python3
"""M7c-9: diagnose Arnoldi acceptance against an independent action oracle.

No production step doubling. The sparse exponential action is diagnostic only.
This is not a certified estimator and does not alter the production solver.
"""
import argparse
import csv
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from scipy.sparse import bmat, csc_matrix, eye
from scipy.sparse.linalg import expm_multiply

from global_error_transport import ROOT, load_operators
from global_error_adaptive_krylov import initial_state

FIELDS = ('n','production_theta','dt_s','probe','requested_rtol','step',
          'action_id','depth','time_s','dimension','accepted','reason',
          'full_input_norm','first_input_norm','first_oracle_norm',
          'full_oracle_norm','indicator','indicator_over_full_scale',
          'indicator_over_first_scale','checkpoint_difference',
          'actual_full_action_error','actual_first_action_error',
          'actual_first_relative_error','jacobian_products_so_far')


def diagnostic_action(a, v, t, rtol, dims, counts, rows, metadata,
                      depth=0, max_depth=16):
    """Replicate M7c-8 decisions, recording true action errors at checkpoints."""
    beta = np.linalg.norm(v)
    size = a.shape[0] // 3
    action_id = counts['actions']
    counts['actions'] += 1
    # Independent reference: sparse expm_multiply, not the Arnoldi controller.
    reference = expm_multiply(t * a, v)
    full_ref = np.linalg.norm(reference)
    first_ref = np.linalg.norm(reference[:size])
    if beta == 0:
        return v.copy()
    cap = min(dims[-1], len(v))
    q = np.zeros((len(v),cap+1))
    hh = np.zeros((cap+1,cap))
    q[:,0] = v / beta
    previous = None
    for k in range(cap):
        w = a @ q[:,k]
        counts['jv'] += 1
        for _ in range(2):
            for i in range(k+1):
                c = np.dot(q[:,i],w)
                hh[i,k] += c
                w -= c*q[:,i]
        hnext = np.linalg.norm(w)
        hh[k+1,k] = hnext
        m = k+1
        breakdown = hnext <= 1e-13 * max(1.,np.linalg.norm(hh[:m,:m],ord=np.inf))
        if not breakdown and m < cap:
            q[:,m] = w / hnext
        if m not in dims and not breakdown and m != cap:
            continue
        small = np.zeros((m+1,m+1))
        small[:m,:m] = hh[:m,:m]
        small[m,m-1] = 0. if breakdown else hnext
        col = beta * expm(t*small)[:,0]
        current = q[:,:m] @ col[:m]
        indicator = abs(col[m])
        scale = max(np.linalg.norm(current),beta,1e-30)
        ratio = indicator / scale
        difference = (np.linalg.norm(current-previous)/scale
                      if previous is not None else np.inf)
        accept = breakdown or (previous is not None and ratio <= rtol
                               and difference <= rtol)
        if breakdown:
            reason = 'breakdown'
        elif accept:
            reason = 'threshold'
        elif m == cap:
            reason = 'subdivide'
        else:
            reason = 'continue'
        err = current-reference
        row = dict(metadata,action_id=action_id,depth=depth,time_s=t,
                   dimension=m,accepted=int(accept),reason=reason,
                   full_input_norm=beta,first_input_norm=np.linalg.norm(v[:size]),
                   first_oracle_norm=first_ref,full_oracle_norm=full_ref,
                   indicator=indicator,indicator_over_full_scale=ratio,
                   indicator_over_first_scale=indicator/max(first_ref,1e-30),
                   checkpoint_difference=difference,
                   actual_full_action_error=np.linalg.norm(err),
                   actual_first_action_error=np.linalg.norm(err[:size]),
                   actual_first_relative_error=np.linalg.norm(err[:size])/max(first_ref,1e-30),
                   jacobian_products_so_far=counts['jv'])
        rows.append(row)
        if accept:
            return current
        previous = current
    if depth >= max_depth:
        raise RuntimeError('Exceeded maximum action subdivision depth')
    mid = diagnostic_action(a,v,t/2,rtol/2,dims,counts,rows,metadata,
                            depth+1,max_depth)
    return diagnostic_action(a,mid,t/2,rtol/2,dims,counts,rows,metadata,
                             depth+1,max_depth)


def run_case(j,n,theta,dt,duration,tol,dims,probe,seed,max_depth,rows):
    size = len(j)
    ident = np.eye(size)
    g = lu_solve(lu_factor(ident-theta*dt*j),ident+(1-theta)*dt*j)
    e = expm(dt*j) # reference only
    js = csc_matrix(j)
    z = csc_matrix((size,size))
    ii = eye(size,format='csc')
    a = bmat([[js,ii,z],[z,z,ii],[z,z,z]],format='csc')
    steps = round(duration/dt)
    if not np.isclose(steps*dt,duration):
        raise ValueError('duration must be an integer multiple of dt')
    y = initial_state(j,n,probe,seed)
    exact = y.copy()
    estimate = np.zeros(size)
    counts = {'actions':0,'jv':0}
    for step in range(steps):
        ynext = g@y
        r0 = (ynext-y)/dt-j@y
        r1 = (ynext-y)/dt-j@ynext
        start = np.concatenate((estimate,r0,(r1-r0)/dt))
        metadata = dict(n=n,production_theta=theta,dt_s=dt,
                        probe=probe,requested_rtol=tol,step=step+1)
        estimate = diagnostic_action(a,start,dt,tol,dims,counts,rows,
                                     metadata,max_depth=max_depth)[:size]
        exact = e@exact
        y = ynext
    actual = y-exact
    abs_error = np.linalg.norm(estimate-actual)
    true_norm = np.linalg.norm(actual)
    accepted = [r for r in rows if r['n']==n and r['production_theta']==theta
                and r['dt_s']==dt and r['probe']==probe
                and r['requested_rtol']==tol and r['accepted']]
    worst = max(accepted,key=lambda r:r['actual_first_relative_error'])
    print(f'{n} theta={theta:g} dt={dt:g} {probe} tol={tol:g} '
          f'global_err={abs_error/max(true_norm,1e-14):.5g} '
          f'|true_err|={true_norm:.5g} abs_diff={abs_error:.5g} '
          f'Jv={counts["jv"]} accepted={len(accepted)} '
          f'worst_action_rel={worst["actual_first_relative_error"]:.5g} '
          f'at step={worst["step"]} m={worst["dimension"]}',flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n',type=int,nargs='+',default=[40])
    p.add_argument('--theta',type=float,nargs='+',default=[.65,1.])
    p.add_argument('--dt',type=float,nargs='+',default=[60.])
    p.add_argument('--probes',nargs='+',default=['fast','slow','impulse'],
                   choices=['impulse','random','slow','fast'])
    p.add_argument('--tol',type=float,nargs='+',default=[1e-2,1e-4])
    p.add_argument('--seed',type=int,default=1729)
    p.add_argument('--max-dim',type=int,default=32)
    p.add_argument('--checkpoint',type=int,default=4)
    p.add_argument('--max-depth',type=int,default=16)
    p.add_argument('--duration',type=float,default=600.)
    p.add_argument('--output',type=Path,default=ROOT/'results-m7c'/'output_diagnostic.csv')
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
                        run_case(operators[n],n,theta,dt,args.duration,tol,
                                 dims,probe,args.seed,args.max_depth,rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print('Wrote',args.output,'with',len(rows),'checkpoint records')

if __name__=='__main__':
    main()
