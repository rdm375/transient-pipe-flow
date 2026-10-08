#!/usr/bin/env python3
"""M7c-17: compare cheap residual transport estimates to one-step truth.

Place beside local_indicator_m7c16.py and boundary_cfl_experiment.py.
Uses artificial endpoint forcing, NOT validated physical boundary conditions.
All estimators use the same Gauss-Legendre residual samples. No step doubling.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.linalg import expm, lu_factor, lu_solve
from local_indicator_m7c16 import ROOT, load_operators, forcing, profile, solve_exact, theta_step


def evaluate(j, b, tau, kind, t, h, theta, y, yn, f0, f1, samples,
             reference, scale, exact_transport):
    x, weights = np.polynomial.legendre.leggauss(samples)
    nodes = (x + 1.0) / 2.0
    weights = weights / 2.0
    dy = (yn - y) / h
    rs = np.array([dy - j @ (y + s * (yn-y)) - b * profile(t+s*h, tau, kind)
                   for s in nodes])
    # Residual moments: h * int_0^1 (1-s)^k r(s) ds.
    moments = [h * np.einsum('i,ij->j', weights * (1-nodes)**k, rs)
               for k in range(3)]
    v0 = moments[0]
    v1 = v0 + j @ moments[1] * h
    v2 = v1 + (j @ (j @ moments[2])) * (h*h/2)
    # Shared theta matrix rational approximation, one extra solve:
    # exp((1-s)hJ) ~ I + ((1-s)/theta) * [(I-theta*hJ)^-1 - I]
    # This matches the first-order expansion and reuses the theta-step factor.
    factor = lu_factor(np.eye(len(y)) - theta*h*j)
    shared = v0 + (lu_solve(factor, moments[1]) - moments[1]) / theta
    # Nodewise backward-Euler resolvent (Q distinct factorizations), for comparison.
    be = np.zeros_like(y)
    for s, w, r in zip(nodes, weights, rs):
        be += h*w*lu_solve(lu_factor(np.eye(len(y)) - (1-s)*h*j), r)
    transported = None
    if exact_transport:
        transported = np.zeros_like(y)
        for s, w, r in zip(nodes, weights, rs):
            transported += h*w*(expm((1-s)*h*j) @ r)
    endpoint = abs(theta-0.5)*h*np.linalg.norm(f1-f0)
    l1 = h*np.dot(weights, np.linalg.norm(rs, axis=1))
    truth = yn - reference
    estimators = {'endpoint': None, 'residual_L1': None,
                  'signed': v0, 'taylor1': v1, 'taylor2': v2,
                  'shared_theta_resolvent': shared, 'nodewise_BE': be}
    if transported is not None:
        estimators['exact_exp_quadrature'] = transported
    actual = np.linalg.norm(truth)
    rows=[]
    for name, v in estimators.items():
        magnitude = endpoint if name == 'endpoint' else l1 if name == 'residual_L1' else np.linalg.norm(v)
        rows.append(dict(estimator=name, actual=actual/scale,
                         estimate=magnitude/scale,
                         ratio=magnitude/actual if actual else np.nan,
                         vector_discrepancy=(np.linalg.norm(v-truth)/actual
                                             if v is not None and actual else np.nan),
                         # Approximate additional work beyond theta step; quadrature
                         # residual J@yi calls are counted separately in `samples`.
                         extra_jv={'endpoint':0,'residual_L1':0,'signed':0,
                                   'taylor1':1,'taylor2':2,
                                   'shared_theta_resolvent':0,'nodewise_BE':0,
                                   'exact_exp_quadrature':0}[name],
                         extra_solves={'shared_theta_resolvent':1,'nodewise_BE':samples}.get(name,0),
                         extra_factorizations={'nodewise_BE':samples}.get(name,0),
                         extra_expm={'exact_exp_quadrature':samples}.get(name,0)))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=40)
    p.add_argument('--theta',type=float,default=.65)
    p.add_argument('--tau',nargs='+',type=float,default=[30,120,600])
    p.add_argument('--profiles',nargs='+',choices=['smooth','rapid','sinusoid'],default=['smooth','sinusoid'])
    p.add_argument('--dt',nargs='+',type=float,default=[7.5,15,30,60])
    p.add_argument('--samples',nargs='+',type=int,default=[8,24])
    p.add_argument('--start-fractions',nargs='+',type=float,default=[0,.25,.5])
    p.add_argument('--duration',type=float,default=600)
    p.add_argument('--exact-transport',action='store_true',help='expensive quadrature oracle')
    p.add_argument('--side',choices=['left','right'],default='left')
    p.add_argument('--component',type=int,choices=[0,1],default=1)
    p.add_argument('--output-dir',type=Path,default=ROOT/'results-m7c')
    a=p.parse_args()
    if a.theta<=0 or min(a.dt)<=0 or min(a.samples)<1 or min(a.tau)<=0 or any(not 0<=v<1 for v in a.start_fractions):
        p.error('theta, dt, samples and tau must be positive; start fractions in [0,1)')
    j=load_operators()[a.n]
    b=forcing(a.n,a.side,a.component)
    rows=[]
    for tau in a.tau:
        for kind in a.profiles:
            ref=solve_exact(j,b,tau,kind,0,a.duration,np.zeros(len(j)),min(tau/40,a.duration/300))
            scale=max(1e-14,max(np.linalg.norm(ref.sol(tt)) for tt in np.linspace(0,a.duration,501)))
            for frac in a.start_fractions:
                t=frac*a.duration
                y=ref.sol(t)
                for h in a.dt:
                    if t+h>a.duration+1e-10: continue
                    yn,f0,f1=theta_step(j,b,tau,kind,t,h,a.theta,y)
                    oracle=solve_exact(j,b,tau,kind,t,t+h,y,min(h/40,tau/40)).y[:,-1]
                    for q in a.samples:
                        result=evaluate(j,b,tau,kind,t,h,a.theta,y,yn,f0,f1,q,oracle,scale,a.exact_transport)
                        for row in result:
                            row.update(n=a.n,theta=a.theta,profile=kind,tau_s=tau,start_s=t,dt_s=h,samples=q)
                            rows.append(row)
                        best=min((r for r in result if np.isfinite(r['vector_discrepancy'])),
                                 key=lambda r:r['vector_discrepancy'])
                        print(f'{kind:8s} tau={tau:5g} t={t:6g} h={h:5g} Q={q:2d} '
                              f'actual={result[0]["actual"]:.3e} best={best["estimator"]} '
                              f'D={best["vector_discrepancy"]:.3e}',flush=True)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    path=a.output_dir/'residual_transport_m7c17.csv'
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print('Wrote',path)
    print('NOTE: shared-theta resolvent counts a solve but assumes existing factorization can be reused.')
    print('NOTE: endpoint and residual_L1 are magnitude-only; vector discrepancy undefined.')
    print('NOTE: residual sampling cost = Q matrix-vector evaluations for all residual-based estimators.')
    print('NOTE: forcing is an artificial endpoint proxy, not validated physical BC.')

if __name__=='__main__': main()
