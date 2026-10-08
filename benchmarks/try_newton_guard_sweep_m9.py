#!/usr/bin/env python3
"""M9 guarded Newton threshold and grid sensitivity, without modifying production sources."""
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
    p.add_argument('--n', type=int, nargs='+', default=[20])
    p.add_argument('--scaled-tol', type=float, nargs='+', default=[1e-13, 1e-12, 1e-11])
    p.add_argument('--dt', type=float, nargs='+', default=[1.875, 3.75, 7.5])
    p.add_argument('--duration', type=float, default=3600)
    p.add_argument('--outdir', type=Path, default=Path('benchmarks/results-m7c/inventory-m7c18a/newton-diagnostic'))
    a = p.parse_args()
    root = Path.cwd()
    original = (root / 'src/fortran77/newton_solver.f').read_text()
    needle = '         INFO = 4\n         RETURN'
    if original.count(needle) != 1:
        p.error('Could not uniquely identify INFO=4 line-search failure branch')
    diagnostics = """         RC = 0D0
         RM = 0D0
         DO 210 I=2,NU,2
            RC = DMAX1(RC,DABS(R(I)))
  210    CONTINUE
         DO 220 I=3,NU-1,2
            RM = DMAX1(RM,DABS(R(I)))
  220    CONTINUE
         DP = 0D0
         DM = DABS(DU(2))
         DO 230 I=1,NU
            IF (I .EQ. 1 .OR. MOD(I,2) .EQ. 0) THEN
               IF (I .NE. 2) DP=DMAX1(DP,DABS(DU(I)))
            ELSE
               DM=DMAX1(DM,DABS(DU(I)))
            ENDIF
  230    CONTINUE
         PSC = DMAX1(1D0,DABS(PINN))
         MSC = DMAX1(1D0,DABS(MOUTN),DABS(MOUTO))
         IF (RNORM .LE. 10D0*RTOL .AND.
     &       RC/MSC .LE. 1D-12 .AND.
     &       RM/(PSC/DX) .LE. 1D-12 .AND.
     &       DABS(R(1))/PSC .LE. 1D-12 .AND.
     &       DP/PSC .LE. 1D-12 .AND.
     &       DM/MSC .LE. 1D-12) THEN
            WRITE(*,*) 'M9 GUARDED STAGNATION ACCEPT',
     &          ' RC=',RC,' RM=',RM,' DP=',DP,' DM=',DM
            NITER = ITER + 1
            INFO = 0
            RETURN
         ENDIF
         WRITE(*,*) 'M9 GUARDED STAGNATION REJECT',
     &       ' RNORM=',RNORM,' RC=',RC,' RM=',RM,
     &       ' DP=',DP,' DM=',DM
         NITER = ITER + 1
         INFO = 4
         RETURN"""
    original = original.replace(
        '      DOUBLE PRECISION RNORM,RNEW,SNORM,USCALE,LAMBDA',
        '      DOUBLE PRECISION RNORM,RNEW,SNORM,USCALE,LAMBDA\n'
        '      DOUBLE PRECISION RC,RM,DP,DM,PSC,MSC',
        1,
    )

    a.outdir.mkdir(parents=True, exist_ok=True)
    results = []
    with tempfile.TemporaryDirectory(prefix='m9_guard_sweep_') as td:
        tmp = Path(td)
        for tol in a.scaled_tol:
            if not (0 < tol < 1):
                p.error('--scaled-tol must be between zero and one')
            tol_literal = f'{tol:.16E}'.replace('E', 'D')
            diagnostic_tol = diagnostics.replace('1D-12', tol_literal)
            instrumented = tmp / 'newton_solver_diagnostic.f'
            instrumented.write_text(original.replace(needle, diagnostic_tol))
            sources = [str(instrumented) if s.endswith('/newton_solver.f')
                       else str(root / s) for s in SOURCES]
            exe = tmp / 'inventory_history_diagnostic'
            cmd = ['gfortran', '-O2', '-ffixed-line-length-none', '-o', str(exe),
                   *sources, str(root / 'benchmarks/inventory_history_m7c18a.f')]
            subprocess.run(cmd, check=True)
            for n in a.n:
                for dt in a.dt:
                    tag = f'n{n}_dt{dt:g}_tol{tol:.0e}'
                    outfile = a.outdir / f'{tag}.csv'
                    proc = subprocess.run(
                        [str(exe), str(dt), str(a.duration), str(outfile), str(n)],
                        capture_output=True, text=True)
                    logtext = proc.stdout + proc.stderr
                    (a.outdir / f'{tag}.log').write_text(logtext)
                    accepts = logtext.count('M9 GUARDED STAGNATION ACCEPT')
                    rejects = logtext.count('M9 GUARDED STAGNATION REJECT')
                    last = ''
                    if outfile.exists():
                        with outfile.open() as f:
                            next(f, None)
                            last = next(reversed(f.readlines()), '').strip()
                    results.append((n, dt, tol, proc.returncode, accepts, rejects, last))
                    print(f'n={n:3d} dt={dt:7g} tol={tol:.0e} '
                          f'rc={proc.returncode} accepts={accepts} rejects={rejects}',
                          flush=True)
    import csv
    with (a.outdir / 'summary.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 'dt', 'scaled_tol', 'returncode', 'guard_accepts',
                    'guard_rejects', 'last_csv_row'])
        w.writerows(results)
    print(f'\nSummary: {a.outdir / "summary.csv"}')
    print('\nProduction solver unchanged. Experimental convergence rule only.')

if __name__ == '__main__':
    main()
