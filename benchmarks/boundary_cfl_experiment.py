#!/usr/bin/env python3
"""M7c-14 exploratory boundary-forcing/CFL study using frozen M7b Jacobians.

IMPORTANT: Endpoint forcing vectors are *proxies*, not derived physical
boundary-condition operators. Supply L and c for an illustrative acoustic CFL.
No step doubling. Independent solve_ivp reference uses DOP853.
"""
import argparse
import csv
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, lu_factor, lu_solve
from global_error_transport import ROOT, load_operators
from global_error_budgeting import candidates

def profile(t, tau, kind):
    if kind == 'smooth':
        x = np.clip(t / tau, 0.0, 1.0)
        return x*x*(3-2*x)
    if kind == 'rapid':
        x = np.clip(t / tau, 0.0, 1.0)
        return x
    if kind == 'sinusoid':
        return np.sin(2*np.pi*t/tau)
    raise ValueError(kind)

def forcing(n, side, component):
    # Interleaved state assumption: [state0_cell0, state1_cell0, ...].
    # The precise physical interpretation of each component MUST be
    # checked against the Fortran state layout before physical conclusions.
    b=np.zeros(2*n)
    b[2*(0 if side=='left' else n-1)+component]=1.
    return b

def simulate(j,n,theta,dt,duration,tau,kind,b,rtol,atol):
    steps=int(round(duration/dt))
    if not np.isclose(steps*dt,duration):
        raise ValueError('duration/dt must be an integer')
    ident=np.eye(len(j))
    lhs=lu_factor(ident-theta*dt*j)
    y=np.zeros(len(j))
    numerical=[y.copy()]
    for k in range(steps):
        t=k*dt
        f0=profile(t,tau,kind)*b
        f1=profile(t+dt,tau,kind)*b
        y=lu_solve(lhs,(ident+(1-theta)*dt*j)@y+
                   dt*((1-theta)*f0+theta*f1))
        numerical.append(y.copy())
    times=np.linspace(0,duration,steps+1)
    reference=solve_ivp(lambda t,z:j@z+profile(t,tau,kind)*b,
                        (0,duration),np.zeros(len(j)),method='DOP853',
                        t_eval=times,rtol=rtol,atol=atol,
                        max_step=min(dt/4,tau/12))
    if not reference.success: raise RuntimeError(reference.message)
    return times,np.asarray(numerical),reference.y.T,reference.nfev

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n',type=int,nargs='+',default=[20,40])
    ap.add_argument('--theta',type=float,default=.65)
    ap.add_argument('--dt',type=float,nargs='+',default=[7.5,15.,30.,60.])
    ap.add_argument('--tau',type=float,nargs='+',default=[30.,120.,600.])
    ap.add_argument('--profiles',nargs='+',default=['smooth','rapid','sinusoid'])
    ap.add_argument('--duration',type=float,default=600.)
    ap.add_argument('--length-m',type=float,required=True,
                    help='Physical pipe length; must match the matrix model')
    ap.add_argument('--sound-speed-ms',type=float,required=True,
                    help='Isothermal sound speed; must match matrix model')
    ap.add_argument('--flow-speed-ms',type=float,default=0.,
                    help='Magnitude of characteristic base flow velocity')
    ap.add_argument('--side',choices=['left','right'],default='left')
    ap.add_argument('--component',type=int,choices=[0,1],default=1)
    ap.add_argument('--amplitude',type=float,default=1.,
                    help='Magnitude of forcing vector (not boundary pressure units)')
    ap.add_argument('--krylov-dims',type=int,nargs='+',default=[4,8,12,16,20,24,28,32])
    ap.add_argument('--rtol',type=float,default=1e-10)
    ap.add_argument('--atol',type=float,default=1e-12)
    ap.add_argument('--output-dir',type=Path,default=ROOT/'results-m7c')
    args=ap.parse_args()
    if args.length_m<=0 or args.sound_speed_ms<=0: ap.error('L and c must be positive')
    operators=load_operators()
    summary=[]; modes=[]
    for n in args.n:
        j=operators[n]
        b=forcing(n,args.side,args.component)*args.amplitude
        dx=args.length_m/n
        eigenvalues,eigenvectors=np.linalg.eig(j)
        # Complex modal decomposition may be ill-conditioned for nonnormal J.
        modal_condition=np.linalg.cond(eigenvectors)
        for tau in args.tau:
            for kind in args.profiles:
                for dt in args.dt:
                    start=time.perf_counter()
                    times,num,ref,nfev=simulate(j,n,args.theta,dt,args.duration,
                                                tau,kind,b,args.rtol,args.atol)
                    difference=num-ref
                    abs_err=np.linalg.norm(difference[-1])
                    scale=max(np.linalg.norm(ref[-1]),1e-14)
                    peak_err=max(np.linalg.norm(v) for v in difference)
                    peak_scale=max(max(np.linalg.norm(v) for v in ref),1e-14)
                    # A physically driven state, not an arbitrary impulse:
                    sample=ref[len(ref)//2]
                    exp_reference=expm(dt*j)@sample
                    krylov=[]
                    if np.linalg.norm(sample)>1e-14:
                        results=candidates(j,sample,dt,args.krylov_dims,len(j))
                        for m in args.krylov_dims:
                            val=results[m]
                            rel=np.linalg.norm(val-exp_reference)/max(
                                np.linalg.norm(exp_reference),1e-14)
                            krylov.append((m,rel))
                    cfl=(args.sound_speed_ms+abs(args.flow_speed_ms))*dt/dx
                    record=dict(n=n,theta=args.theta,dt_s=dt,tau_s=tau,
                                profile=kind,cfl=cfl,forcing_ratio=dt/tau,
                                acoustic_to_forcing=args.length_m/args.sound_speed_ms/tau,
                                final_rel_error=abs_err/scale,
                                peak_rel_error=peak_err/peak_scale,
                                reference_nfev=nfev,modal_cond=modal_condition,
                                first_krylov_dim_1pct=next((m for m,v in krylov if v<=.01),''),
                                first_krylov_dim_1e4=next((m for m,v in krylov if v<=1e-4),''),
                                elapsed_ms=1000*(time.perf_counter()-start))
                    summary.append(record)
                    print(f'n={n:2d} {kind:8s} tau={tau:6g}s dt={dt:5g}s '
                          f'CFL={cfl:7.3f} final={record["final_rel_error"]:.3e} '
                          f'peak={record["peak_rel_error"]:.3e} '
                          f'm_1pct={record["first_krylov_dim_1pct"]} '
                          f'm_1e4={record["first_krylov_dim_1e4"]}',flush=True)
                    for m,v in krylov:
                        modes.append(dict(n=n,dt_s=dt,tau_s=tau,profile=kind,
                                          cfl=cfl,krylov_dim=m,action_rel_error=v))
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for name,rows in [('boundary_cfl_summary.csv',summary),
                      ('boundary_cfl_krylov.csv',modes)]:
        path=args.output_dir/name
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
        print('Wrote',path)
    print('CAUTION: endpoint forcing is a state-space proxy, not a verified physical BC.')
    print('CFL is based on supplied L/c; verify those match M7b matrix provenance.')
    print('Krylov errors are for ONE frozen-state exponential action, not global-error estimation.')
if __name__=='__main__':
    main()
