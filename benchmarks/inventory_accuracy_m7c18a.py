#!/usr/bin/env python3
"""M7c-18A: Fortran inventory accuracy versus timestep; no step doubling."""
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile
import numpy as np

SOURCES = [
    'src/fortran77/models/eos_constant_z.f',
    'src/fortran77/models/friction_swamee_jain.f',
    'src/fortran77/transient_residual.f',
    'src/fortran77/transient_jacobian.f',
    'src/fortran77/linear_solve.f',
    'src/fortran77/newton_solver.f',
    'src/fortran77/transient_step.f',
    'src/fortran77/transient_integrate.f',
]

def load(path):
    with path.open(newline='') as f:
        return {k: np.array([float(r[k]) for r in csv.DictReader(f)]) for k in []}

def read(path):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f))
    return {key: np.asarray([float(row[key]) for row in rows]) for key in rows[0]}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, nargs='+', default=[20, 40])
    ap.add_argument('--dt', type=float, nargs='+', default=[7.5, 15, 30, 60, 120])
    ap.add_argument('--reference-dt', type=float, default=1.875)
    ap.add_argument('--duration', type=float, default=3600)
    ap.add_argument('--outdir', type=Path, default=Path('benchmarks/results-m7c/inventory-m7c18a'))
    args = ap.parse_args()
    root = Path.cwd()
    source = root / 'benchmarks' / 'inventory_history_m7c18a.f'
    missing = [str(p) for p in [source] + [root / s for s in SOURCES] if not p.exists()]
    if missing:
        ap.error('Missing files: ' + ', '.join(missing))
    args.outdir.mkdir(parents=True, exist_ok=True)
    exe = args.outdir / 'inventory_history_m7c18a'
    subprocess.run(['gfortran', '-O2', '-ffixed-line-length-none', '-o', str(exe),
                    *SOURCES, str(source)], check=True)
    columns = ['n', 'dt_s', 'reference_dt_s', 'duration_s', 'steps', 'acoustic_cfl',
               'max_inventory_error_kg', 'max_inventory_change_error_kg',
               'final_inventory_change_error_kg', 'max_imbalance_error_kg_s',
               'max_inlet_flow_error_kg_s', 'max_outlet_pressure_error_pa',
               'max_step_conservation_defect_kg', 'cumulative_conservation_defect_kg',
               'max_inventory_change_reference_kg', 'max_reference_imbalance_kg_s']
    results = []
    for n in args.n:
        def run(dt):
            steps = args.duration / dt
            if not math.isclose(steps, round(steps), abs_tol=1e-8):
                raise ValueError(f'duration {args.duration} not divisible by dt {dt}')
            path = args.outdir / f'n{n}_dt{dt:g}.csv'
            subprocess.run([str(exe), str(dt), str(args.duration), str(path), str(n)], check=True)
            return read(path)
        ref = run(args.reference_dt)
        for dt in args.dt:
            cur = run(dt)
            t = cur['time_s']
            def diff(key):
                return cur[key] - np.interp(t, ref['time_s'], ref[key])
            inv = diff('inventory_kg')
            # Both runs start at the same discrete steady profile.
            inv_change = inv - inv[0]
            ref_change = np.interp(t, ref['time_s'], ref['inventory_kg']) - ref['inventory_kg'][0]
            row = dict(
                n=n, dt_s=dt, reference_dt_s=args.reference_dt, duration_s=args.duration,
                steps=len(t)-1,
                acoustic_cfl=dt * math.sqrt(0.9 * 500.0 * 288.15) / (100000.0 / n),
                max_inventory_error_kg=np.max(np.abs(inv)),
                max_inventory_change_error_kg=np.max(np.abs(inv_change)),
                final_inventory_change_error_kg=inv_change[-1],
                max_imbalance_error_kg_s=np.max(np.abs(diff('imbalance_kg_s'))),
                max_inlet_flow_error_kg_s=np.max(np.abs(diff('inlet_kg_s'))),
                max_outlet_pressure_error_pa=np.max(np.abs(diff('outlet_pressure_pa'))),
                max_step_conservation_defect_kg=np.max(np.abs(cur['step_balance_defect_kg'][1:])),
                cumulative_conservation_defect_kg=cur['cumulative_defect_kg'][-1],
                max_inventory_change_reference_kg=np.max(np.abs(ref_change)),
                max_reference_imbalance_kg_s=np.max(np.abs(np.interp(t, ref['time_s'], ref['imbalance_kg_s']))),
            )
            results.append(row)
            print(f'n={n:3d} dt={dt:7g} CFL={row["acoustic_cfl"]:6.2f} '
                  f'max |dM error|={row["max_inventory_change_error_kg"]:11.4g} kg '
                  f'max |imbalance error|={row["max_imbalance_error_kg_s"]:10.4g} kg/s '
                  f'max |balance defect|={row["max_step_conservation_defect_kg"]:.3g} kg', flush=True)
    with (args.outdir / 'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(results)
    print('Wrote', args.outdir / 'summary.csv')
    print('Reference is a finer theta-method run, not an exact solution; check convergence by halving --reference-dt.')
    print('Mass conservation and inventory accuracy are distinct diagnostics.')

if __name__ == '__main__':
    main()
