#!/usr/bin/env python3
"""Compare nonlinear F77 responses with production-Jacobian linear theta model.

The forcing is the prescribed outlet mass flow. For fixed inlet pressure,
its tangent forcing enters only the outlet half-cell continuity equation.
"""
import csv
import math
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / 'results-m7b'
AMPS = (1.0, 0.5, 0.25, 0.125)
THETAS = (0.5, 0.65, 1.0)
DTS = (15.0, 60.0)
N = 40


def load_operator():
    mat = np.zeros((2*N, 2*N))
    masses = np.zeros(2*N)
    with (OUT/'modal_matrix.csv').open(newline='') as f:
        for r in csv.DictReader(f):
            if int(r['n']) != N:
                continue
            i, j = int(r['row'])-1, int(r['col'])-1
            mat[i, j] = float(r['operator'])
            masses[i] = float(r['mass'])
    assert np.all(masses > 0)
    return mat, masses


def demand(t):
    s = min(1.0, max(0.0, t/600.0))
    return 0.5*(1.0-math.cos(math.pi*s))


def linear_history(mat, masses, dt, theta):
    # Semidiscrete tangent: y'=A y + b q(t).
    # Last equation is outlet half-cell continuity, so b=-1/M_last.
    b = np.zeros(2*N)
    b[-1] = -1.0/masses[-1]
    eye = np.eye(2*N)
    left = eye-theta*dt*mat
    prop = np.linalg.solve(left, eye+(1.0-theta)*dt*mat)
    drive = np.linalg.solve(left, dt*b)
    y = np.zeros(2*N)
    rows = {}
    for k in range(int(round(1800/dt))+1):
        t = k*dt
        # Four nodal pressure sensors, inlet flow, and total linepack.
        # The inlet half-cell pressure perturbation is zero.
        dx = 100000.0/N
        area = math.pi/4
        crho = 1.0/(0.90*500.0*288.15)
        weights = np.full(N, area*dx*crho)
        weights[-1] *= 0.5
        linepack = np.dot(weights, y[1::2])
        rows[round(t, 9)] = np.array([
            y[2*(N//4)-1], y[2*(N//2)-1],
            y[2*(3*N//4)-1], y[-1], y[0], linepack
        ])
        if k < int(round(1800/dt)):
            y = prop@y + drive*(theta*demand(t+dt)
                                  +(1-theta)*demand(t))
    return rows


def main():
    mat, masses = load_operator()
    assert np.all(np.isfinite(mat))
    reference = {}
    for theta in THETAS:
        for dt in DTS:
            reference[(dt,theta)] = linear_history(mat,masses,dt,theta)
    fields = ('p25_pa','p50_pa','p75_pa','pout_pa','min_kg_s','linepack_kg')
    raw = {}
    with (OUT/'small_nonlinear.csv').open(newline='') as f:
        for r in csv.DictReader(f):
            case = int(r['case_id'])
            raw.setdefault(case, []).append(r)
    assert len(raw)==24, f'expected 24 histories, got {len(raw)}'
    results=[]
    for ia,amp in enumerate(AMPS):
        for it,theta in enumerate(THETAS):
            for ih,dt in enumerate(DTS):
                case = ia*6+it*2+ih+1
                history=raw[case]
                base=np.array([float(history[0][key]) for key in fields])
                error=[]
                # Compare perturbations, not absolute state values.
                # Pressure and linepack are scaled into their natural response
                # magnitudes; retain raw errors for dimensional diagnostics.
                for r in history:
                    t=round(float(r['time']),9)
                    actual=np.array([float(r[key]) for key in fields])-base
                    predicted=amp*reference[(dt,theta)][t]
                    error.append(actual-predicted)
                errors=np.array(error)
                max_p=np.max(np.abs(errors[:,:4]))
                max_m=np.max(np.abs(errors[:,4]))
                max_l=np.max(np.abs(errors[:,5]))
                results.append(dict(amplitude_kg_s=amp,dt_s=dt,theta=theta,
                    max_pressure_discrepancy_pa=max_p,
                    max_inlet_flow_discrepancy_kg_s=max_m,
                    max_linepack_discrepancy_kg=max_l))
    with (OUT/'small_metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]))
        writer.writeheader();writer.writerows(results)
    report=['# M7b — Small-perturbation linearization verification','',
      'The linear model uses the production analytical Jacobian and the',
      'outlet half-cell forcing derivative. Both models use the same',
      'theta-method and exactly aligned smooth demand forcing.','',
      'The comparison is made at matching timesteps, so this test examines',
      'linearization consistency, **not temporal discretization accuracy**.','',
      '| dt (s) | theta | pressure error at amplitude 1 (Pa) | pressure error at amplitude 0.125 (Pa) | observed order |',
      '|---:|---:|---:|---:|---:|']
    orders=[]
    for dt in DTS:
        for theta in THETAS:
            sub=[r for r in results if r['dt_s']==dt and r['theta']==theta]
            e0=sub[0]['max_pressure_discrepancy_pa']
            e1=sub[-1]['max_pressure_discrepancy_pa']
            order=math.log(e0/e1)/math.log(8) if e1>0 else math.inf
            orders.append(order)
            report.append(f'| {dt:g} | {theta:.2f} | {e0:.7g} | {e1:.7g} | {order:.3f} |')
    report += ['', '## Interpretation', '',
      '- Quadratic discrepancy in perturbation amplitude is expected from',
      '  a smooth nonlinear model with a consistent first-order tangent.',
      '- Very small perturbations eventually encounter roundoff and Newton',
      '  convergence tolerance floors; measured orders need not be exactly 2.',
      '- Agreement verifies the linearized forced time-step response, not',
      '  a decomposition of finite-amplitude histories into individual modes.',
      '- The eigenvalues of the same operator independently determine the',
      '  theta amplification factors reported in `MODAL_REPORT.md`.',
      '- This does not verify continuum wave speed or physical validity.','']
    (OUT/'SMALL_REPORT.md').write_text('\n'.join(report))
    print('Observed pressure amplitude orders:', [round(x,3) for x in orders])
    if min(orders)<1.5:
        raise SystemExit('FAIL: small-perturbation pressure error not near quadratic')
    print('PASS: nonlinear/linear small-perturbation consistency')

if __name__=='__main__':
    main()
