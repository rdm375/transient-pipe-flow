#!/usr/bin/env python3
"""M7c-15: fixed/boundary/defect/hybrid timestep comparison.

Frozen linearized M7b Jacobian, artificial endpoint state forcing (NOT a
physical BC operator). No step doubling. Defect is the theta-method's
endpoint derivative mismatch, O(h^2) for theta != 1/2. Reference DOP853.
"""
import argparse
import csv
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import lu_factor, lu_solve
from global_error_transport import ROOT, load_operators
from boundary_cfl_experiment import profile, forcing


def reference_solution(j, b, tau, kind, duration, rtol=1e-11, atol=1e-13):
    # Avoid solver steps skipping oscillations, and resolve forcing corners.
    max_step = min(tau/30, duration/200)
    sol = solve_ivp(lambda t, y: j @ y + b * profile(t, tau, kind),
                    (0, duration), np.zeros(len(j)), method='DOP853',
                    rtol=rtol, atol=atol, max_step=max_step, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def solve_case(j, b, tau, kind, duration, theta, policy, fixed_h, hmax,
               boundary_fraction, local_tol, reference, reference_scale,
               wave_hmax, min_h=1e-6):
    ident = np.eye(len(j))
    t = 0.0
    y = np.zeros(len(j))
    h = fixed_h if policy == 'fixed' else hmax
    accepts = rejects = factorizations = 0
    max_err = 0.0
    max_cfl = 0.0
    h_min = float('inf')
    h_max = 0.0
    defect_max = 0.0
    start = time.perf_counter()
    while t < duration - 1e-10:
        if accepts + rejects > 100000:
            raise RuntimeError('step limit exceeded')
        h = min(h, duration-t, hmax)
        if policy == 'fixed':
            h = min(fixed_h, duration-t)
        if policy in ('boundary', 'hybrid'):
            # A priori forcing-bandwidth cap. For a ramp, limit during the
            # ramp and land exactly on its terminal kink at t=tau.
            if kind == 'sinusoid':
                h = min(h, boundary_fraction*tau)
            elif t < tau - 1e-10:
                h = min(h, boundary_fraction*tau, tau-t)
        if policy == 'hybrid' and wave_hmax is not None:
            h = min(h, wave_hmax)
        if h < min_h:
            raise RuntimeError(f'h={h} below min_h at t={t}')
        f0 = j@y + b*profile(t,tau,kind)
        f1b = b*profile(t+h,tau,kind)
        lhs = lu_factor(ident-theta*h*j)
        factorizations += 1
        yn = lu_solve(lhs, y+h*((1-theta)*f0+theta*f1b))
        f1 = j@yn+f1b
        # Leading theta quadrature defect (exact-endpoint sign irrelevant).
        # For linear systems f1-f0 is available after the single theta solve.
        defect = abs(theta-.5)*h*np.linalg.norm(f1-f0)/reference_scale
        if policy in ('defect', 'hybrid') and defect > local_tol:
            rejects += 1
            h *= max(.15, .85*np.sqrt(local_tol/max(defect,1e-300)))
            continue
        t += h
        y = yn
        accepts += 1
        max_err = max(max_err,np.linalg.norm(y-reference.sol(t))/reference_scale)
        defect_max = max(defect_max,defect)
        h_min = min(h_min,h)
        h_max = max(h_max,h)
        if policy in ('defect','hybrid'):
            h = min(hmax, h*min(2.,max(.5,.9*np.sqrt(local_tol/max(defect,1e-300)))))
        elif policy == 'boundary':
            h = hmax
    return dict(steps=accepts,rejected=rejects,linear_solves=factorizations,
                max_rel_error=max_err,
                final_rel_error=np.linalg.norm(y-reference.sol(duration))/reference_scale,
                min_dt=h_min,max_dt=h_max,max_local_indicator=defect_max,
                elapsed_ms=(time.perf_counter()-start)*1000)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n',type=int,nargs='+',default=[40])
    ap.add_argument('--theta',type=float,default=.65)
    ap.add_argument('--tau',type=float,nargs='+',default=[30,120,600])
    ap.add_argument('--profiles',nargs='+',choices=['smooth','rapid','sinusoid'],
                    default=['smooth','sinusoid'])
    ap.add_argument('--duration',type=float,default=600)
    ap.add_argument('--fixed-dt',type=float,nargs='+',default=[7.5,15,30,60])
    ap.add_argument('--hmax',type=float,default=60)
    ap.add_argument('--boundary-fraction',type=float,default=.125)
    ap.add_argument('--local-tol',type=float,nargs='+',default=[.001,.005,.02])
    ap.add_argument('--cfl-cap',type=float,default=None,
                    help='Optional acoustic accuracy cap for hybrid; omit to disable')
    ap.add_argument('--length-m',type=float,required=True)
    ap.add_argument('--sound-speed-ms',type=float,required=True)
    ap.add_argument('--flow-speed-ms',type=float,default=0)
    ap.add_argument('--side',choices=['left','right'],default='left')
    ap.add_argument('--component',type=int,choices=[0,1],default=1)
    ap.add_argument('--output-dir',type=Path,default=ROOT/'results-m7c')
    args=ap.parse_args()
    if args.theta == .5:
        ap.error('This O(h^2) indicator vanishes at theta=0.5; choose theta != 0.5')
    if min(args.length_m,args.sound_speed_ms,args.hmax,args.boundary_fraction)<=0:
        ap.error('positive length, speed, hmax, and boundary-fraction required')
    if args.cfl_cap is not None and args.cfl_cap<=0:
        ap.error('positive CFL cap required')
    operators=load_operators()
    rows=[]
    for n in args.n:
        j=operators[n]
        b=forcing(n,args.side,args.component)
        dx=args.length_m/n
        wave_hmax=(args.cfl_cap*dx/(args.sound_speed_ms+abs(args.flow_speed_ms))
                   if args.cfl_cap is not None else None)
        for tau in args.tau:
            for kind in args.profiles:
                ref=reference_solution(j,b,tau,kind,args.duration)
                # Common reference scale, preventing near-zero endpoint denominators.
                check_t=np.linspace(0,args.duration,401)
                scale=max(np.linalg.norm(ref.sol(t)) for t in check_t)
                scale=max(scale,1e-14)
                configs=[('fixed',h,None) for h in args.fixed_dt]
                configs += [('boundary',None,None)]
                for tol in args.local_tol:
                    configs += [('defect',None,tol),('hybrid',None,tol)]
                for policy,fixed_h,tol in configs:
                    result=solve_case(j,b,tau,kind,args.duration,args.theta,
                                      policy,fixed_h,args.hmax,args.boundary_fraction,
                                      tol,ref,scale,wave_hmax)
                    row=dict(n=n,profile=kind,tau_s=tau,policy=policy,
                             requested_fixed_dt=fixed_h if fixed_h is not None else '',
                             local_tol=tol if tol is not None else '',
                             cfl_cap=args.cfl_cap if args.cfl_cap is not None else '',
                             max_actual_cfl=result['max_dt']*(args.sound_speed_ms+
                                              abs(args.flow_speed_ms))/dx,
                             **result)
                    rows.append(row)
                    print(f'n={n} {kind:8s} tau={tau:5g} {policy:8s} '
                          f'h={str(fixed_h or "adaptive"):>8s} '
                          f'tol={str(tol or "-"):>6s} steps={result["steps"]:4d} '
                          f'rej={result["rejected"]:3d} '
                          f'peak={result["max_rel_error"]:.3e} '
                          f'final={result["final_rel_error"]:.3e}',flush=True)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    path=args.output_dir/'adaptive_timestep_m7c15.csv'
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print('Wrote',path)
    print('NOTE: forcing is an endpoint state proxy; not verified physical BC.')
    print('NOTE: errors sampled at accepted endpoints, not continuous-time maxima.')
    print('NOTE: indicator is a leading-order local defect proxy, not an error bound.')

if __name__=='__main__':
    main()
