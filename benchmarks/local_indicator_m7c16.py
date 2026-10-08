#!/usr/bin/env python3
"""M7c-16: audit local theta-step indicators against independent one-step truth.

Uses frozen M7b Jacobians and artificial endpoint forcing (NOT physical BC).
No step doubling. Each candidate starts at an independently integrated exact state.
Compare endpoint defect, sampled residual L1, and propagated residual estimate.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import lu_factor, lu_solve, expm
from global_error_transport import ROOT, load_operators
from boundary_cfl_experiment import profile, forcing


def solve_exact(j, b, tau, kind, t0, t1, y0, max_step):
    out = solve_ivp(lambda t,y: j@y + b*profile(t,tau,kind),
                    (t0,t1), y0, method='DOP853', rtol=2e-11, atol=2e-13,
                    max_step=max_step, dense_output=True)
    if not out.success:
        raise RuntimeError(out.message)
    return out


def theta_step(j,b,tau,kind,t,h,theta,y):
    f0=j@y+b*profile(t,tau,kind)
    rhs=y+h*((1-theta)*f0+theta*b*profile(t+h,tau,kind))
    yn=lu_solve(lu_factor(np.eye(len(j))-theta*h*j),rhs)
    f1=j@yn+b*profile(t+h,tau,kind)
    return yn,f0,f1


def indicators(j,b,tau,kind,t,h,theta,y,yn,f0,f1,samples,propagate):
    dy=(yn-y)/h
    # A sample-based L1 residual can itself alias: deliberately test resolution.
    x,w=np.polynomial.legendre.leggauss(samples)
    s=(x+1)/2; w=w/2
    rs=[]
    for si in s:
        yi=y+si*(yn-y)
        ri=dy-j@yi-b*profile(t+si*h,tau,kind)
        rs.append(ri)
    residual_l1=h*sum(wi*np.linalg.norm(ri) for wi,ri in zip(w,rs))
    # Exact propagation of the interpolant residual, with numerical quadrature.
    # Numerical endpoint minus exact endpoint = integral exp((h-s)J) r(s) ds.
    propagated=np.nan
    if propagate:
        accum=np.zeros_like(y)
        for si,wi,ri in zip(s,w,rs):
            accum+=wi*(expm((1-si)*h*j)@ri)
        propagated=np.linalg.norm(h*accum)
    endpoint=abs(theta-.5)*h*np.linalg.norm(f1-f0)
    return endpoint,residual_l1,propagated


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=40)
    p.add_argument('--theta',type=float,default=.65)
    p.add_argument('--tau',nargs='+',type=float,default=[30,120,600])
    p.add_argument('--profiles',nargs='+',choices=['smooth','rapid','sinusoid'],default=['smooth','sinusoid'])
    p.add_argument('--dt',nargs='+',type=float,default=[7.5,15,30,60])
    p.add_argument('--samples',nargs='+',type=int,default=[3,8,24])
    p.add_argument('--start-fractions',nargs='+',type=float,default=[0.,.25,.5])
    p.add_argument('--duration',type=float,default=600)
    p.add_argument('--propagate',action='store_true',help='Expensive matrix-exponential residual transport')
    p.add_argument('--side',choices=['left','right'],default='left')
    p.add_argument('--component',type=int,choices=[0,1],default=1)
    p.add_argument('--output-dir',type=Path,default=ROOT/'results-m7c')
    a=p.parse_args()
    if min(a.dt)<=0 or min(a.samples)<1 or not 0<=min(a.start_fractions)<=max(a.start_fractions)<1:
        p.error('Require positive dt and samples, start fractions in [0,1)')
    j=load_operators()[a.n]
    b=forcing(a.n,a.side,a.component)
    rows=[]
    for tau in a.tau:
        for kind in a.profiles:
            ref=solve_exact(j,b,tau,kind,0,a.duration,np.zeros(len(j)),min(tau/40,a.duration/300))
            scale=max(np.linalg.norm(ref.sol(tt)) for tt in np.linspace(0,a.duration,501))
            scale=max(scale,1e-14)
            for frac in a.start_fractions:
                t=frac*a.duration
                y=ref.sol(t)
                for h in a.dt:
                    if t+h>a.duration+1e-10:
                        continue
                    yn,f0,f1=theta_step(j,b,tau,kind,t,h,a.theta,y)
                    # Independent one-step oracle starting from exact state at t.
                    truth=solve_exact(j,b,tau,kind,t,t+h,y,min(h/40,tau/40))
                    actual=np.linalg.norm(yn-truth.y[:,-1])/scale
                    for ns in a.samples:
                        endpoint,resid,prop=indicators(j,b,tau,kind,t,h,a.theta,y,yn,f0,f1,ns,a.propagate)
                        row=dict(n=a.n,profile=kind,tau_s=tau,start_s=t,dt_s=h,
                                 samples=ns,theta=a.theta,actual_local=actual,
                                 endpoint_indicator=endpoint/scale,
                                 residual_l1=resid/scale,
                                 propagated_residual=prop/scale if a.propagate else '',
                                 endpoint_ratio=endpoint/(actual*scale) if actual>0 else '',
                                 residual_ratio=resid/(actual*scale) if actual>0 else '')
                        rows.append(row)
                        print(f'{kind:8s} tau={tau:5g} t={t:6g} h={h:5g} '
                              f'Q={ns:2d} actual={actual:.3e} endpoint={endpoint/scale:.3e} '
                              f'L1={resid/scale:.3e}'+
                              (f' transported={prop/scale:.3e}' if a.propagate else ''),flush=True)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    out=a.output_dir/'local_indicator_m7c16.csv'
    with out.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print('Wrote',out)
    print('NOTE: sampled residual quadrature can alias; compare sample counts.')
    print('NOTE: forced endpoint vector is not a verified physical boundary operator.')

if __name__=='__main__':
    main()
