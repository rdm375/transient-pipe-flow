#!/usr/bin/env python3
"""Linearized eigenmode verification from the production F77 Jacobian.

Requires NumPy. This analyzes the constrained semidiscrete operator,
not a continuum PDE dispersion relation or nonlinear modal decomposition.
"""
import csv
import cmath
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent / 'results-m7b'
THETAS = (0.5, 0.65, 1.0)
STEPS = (15., 60., 120.)

def read_matrices():
    matrices = {}
    with (ROOT / 'modal_matrix.csv').open(newline='') as handle:
        for row in csv.DictReader(handle):
            n = int(row['n'])
            if n not in matrices:
                matrices[n] = np.zeros((2*n, 2*n))
            matrices[n][int(row['row'])-1, int(row['col'])-1] = float(row['operator'])
    return matrices

def amplification(z, theta):
    return (1+(1-theta)*z)/(1-theta*z)

def main():
    matrices = read_matrices()
    assert set(matrices) == {20, 40}
    rows=[]
    for n, operator in sorted(matrices.items()):
        assert np.all(np.isfinite(operator))
        eigenvalues=np.linalg.eigvals(operator)
        assert max(x.real for x in eigenvalues) < 1e-8, 'unstable continuous mode'
        for dt in STEPS:
            for theta in THETAS:
                for eig in eigenvalues:
                    g=amplification(dt*eig,theta)
                    assert abs(g) <= 1+1e-10, 'theta amplification unstable'
                    decay=-math.log(abs(g))/dt if abs(g)>0 else math.inf
                    frequency=cmath.phase(g)/dt
                    rows.append(dict(n=n,dt_s=dt,theta=theta,
                        lambda_real_s_inv=eig.real,lambda_imag_rad_s=eig.imag,
                        amplification_abs=abs(g),amplification_phase_rad=cmath.phase(g),
                        discrete_decay_s_inv=decay,
                        continuous_decay_s_inv=-eig.real,
                        discrete_frequency_rad_s=frequency,
                        continuous_frequency_rad_s=eig.imag))
    with (ROOT/'modal_metrics.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    report=['# M7b — Linearized modal verification','',
        'The operator is constructed directly from the production F77 analytical',
        'Jacobian evaluated at the analytical steady state, with fixed inlet',
        'pressure and prescribed outlet flow. The inlet-flow algebraic constraint',
        'is eliminated. The remaining 2N variables are face flows and interior',
        'nodal pressures. Thus the modal results use the **actual discrete**',
        'staggered geometry, including the outlet half-cell.','',
        'For each continuous eigenvalue λ of the semidiscrete operator,',
        'the theta integrator predicts g=(1+(1−θ)hλ)/(1−θhλ).',
        'The decay rate is −log|g|/h and the sampled angular frequency',
        'is arg(g)/h. Frequencies are aliased modulo 2π/h.','',
        '| N | dt (s) | θ | max |g| | min |g| |',
        '|---:|---:|---:|---:|---:|']
    for n in sorted(matrices):
        for dt in STEPS:
            for theta in THETAS:
                subset=[r for r in rows if r['n']==n and r['dt_s']==dt and r['theta']==theta]
                report.append(f"| {n} | {dt:g} | {theta:.2f} | {max(r['amplification_abs'] for r in subset):.6f} | {min(r['amplification_abs'] for r in subset):.6f} |")
    report += ['', '## Interpretation and boundaries','',
        '- Eigenmode damping is a **linearized prediction**, not a fitted nonlinear response.',
        '- The sensor histories include nonlinear friction, finite disturbances, and boundary forcing.',
        '- A slope reversal does not identify a particular eigenmode or prove instability.',
        '- Step/pulse endpoint quadrature can introduce forcing errors distinct from modal damping.',
        '- This calculation does not establish continuum wave speed, a nonlinear error bound, or physical validation.',
        '- A direct small-perturbation experiment is still required to establish modal agreement with nonlinear histories.',
        '- `modal_metrics.csv` contains every complex eigenmode at each N, dt and theta.','']
    (ROOT/'MODAL_REPORT.md').write_text('\n'.join(report))
    print(f'PASS: {len(rows)} modal amplification records; two grids')

if __name__=='__main__':
    main()
