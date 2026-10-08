#!/usr/bin/env python3
"""M7c-18A Newton line-search diagnostics, without modifying production sources."""
import argparse
from pathlib import Path
import subprocess
import tempfile

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

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--n', type=int, default=20)
    p.add_argument('--dt', type=float, nargs='+', default=[1.875, 3.75, 7.5])
    p.add_argument('--duration', type=float, default=3600)
    p.add_argument('--outdir', type=Path, default=Path('benchmarks/results-m7c/inventory-m7c18a/newton-diagnostic'))
    a = p.parse_args()
    root = Path.cwd()
    original = (root / 'src/fortran77/newton_solver.f').read_text()
    needle = '         INFO = 4\n         RETURN'
    if original.count(needle) != 1:
        p.error('Could not uniquely identify INFO=4 line-search failure branch')
    diagnostics = '''         WRITE(*,*) 'M7C18 NEWTON DIAGNOSTIC',
     &       ' ITER=',ITER,' RNORM=',RNORM,
     &       ' RNEW_LAST=',RNEW,' RTOL=',RTOL,
     &       ' SNORM=',SNORM,' USCALE=',USCALE,
     &       ' LAMBDA_LAST=',LAMBDA*2D0
         INFO = 4
         RETURN'''
    a.outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='m7c18_newton_') as td:
        tmp = Path(td)
        instrumented = tmp / 'newton_solver_diagnostic.f'
        instrumented.write_text(original.replace(needle, diagnostics))
        sources = [str(instrumented) if s.endswith('/newton_solver.f') else str(root / s) for s in SOURCES]
        exe = tmp / 'inventory_history_diagnostic'
        cmd = ['gfortran', '-O2', '-ffixed-line-length-none', '-o', str(exe), *sources,
               str(root / 'benchmarks/inventory_history_m7c18a.f')]
        subprocess.run(cmd, check=True)
        for dt in a.dt:
            outfile = a.outdir / f'n{a.n}_dt{dt:g}.csv'
            proc = subprocess.run([str(exe), str(dt), str(a.duration), str(outfile), str(a.n)],
                                  capture_output=True, text=True)
            log = a.outdir / f'n{a.n}_dt{dt:g}.log'
            log.write_text(proc.stdout + proc.stderr)
            print(f'\n=== n={a.n} dt={dt:g} returncode={proc.returncode} ===', flush=True)
            print((proc.stdout + proc.stderr).strip()[-3000:], flush=True)
            print(f'Log: {log}', flush=True)
            if outfile.exists():
                with outfile.open() as f:
                    header = next(f, '')
                    tail = list(f)[-2:]
                print('Last completed timestep(s):')
                for line in tail:
                    print(line.strip(), flush=True)
    print('\nProduction solver unchanged. No tolerances modified.')

if __name__ == '__main__':
    main()
