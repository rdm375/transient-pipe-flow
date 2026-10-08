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
    diagnostics = '''C        M9 continuity backward-error diagnostics.
         WRITE(*,*) 'M9 CONTINUITY BACKWARD ERROR'
         WRITE(*,*) 'ITER=',ITER,' RNORM=',RNORM,
     &              ' RTOL=',RTOL
         WRITE(*,*) 'ROW STORAGE NETFLOW RESIDUAL',
     &              ' ETA_NET ETA_TERMS'

         DO 210 I=0,N
            IF (I .EQ. 0) THEN
               J = 2
               STERM = A*DX/(2D0*DT)*
     &          (EOS_CZ_RHO(U(1),T,Z,RS)
     &          -EOS_CZ_RHO(PO(0),T,Z,RS))
               FIN = THETA*U(2)+(1D0-THETA)*MINO
               FOUT = THETA*U(3)+(1D0-THETA)*MO(0)
            ELSE
               J = 2*I+2
               IF (I .EQ. N) THEN
                  STERM = A*DX/(2D0*DT)*
     &             (EOS_CZ_RHO(U(J),T,Z,RS)
     &             -EOS_CZ_RHO(PO(I),T,Z,RS))
                  FIN = THETA*U(J-1)
     &                 +(1D0-THETA)*MO(I-1)
                  FOUT = THETA*MOUTN
     &                 +(1D0-THETA)*MOUTO
               ELSE
                  STERM = A*DX/DT*
     &             (EOS_CZ_RHO(U(J),T,Z,RS)
     &             -EOS_CZ_RHO(PO(I),T,Z,RS))
                  FIN = THETA*U(J-1)
     &                 +(1D0-THETA)*MO(I-1)
                  FOUT = THETA*U(J+1)
     &                 +(1D0-THETA)*MO(I)
               ENDIF
            ENDIF
            QTERM = FOUT-FIN
            DENNET = DABS(STERM)+DABS(QTERM)
            DENTER = DABS(STERM)+DABS(FIN)+DABS(FOUT)
            ETANET = 0D0
            ETATER = 0D0
            IF (DENNET .GT. 0D0)
     &          ETANET = DABS(R(J))/DENNET
            IF (DENTER .GT. 0D0)
     &          ETATER = DABS(R(J))/DENTER
            WRITE(*,'(I4,5(1X,1PE16.8))')
     &           J,STERM,QTERM,R(J),ETANET,ETATER
  210    CONTINUE

         PMAX = 0D0
         MMAX = 0D0
         PULP = 0D0
         MULP = 0D0
         DO 220 I=1,NU
            IF (I .EQ. 1 .OR. MOD(I,2) .EQ. 0) THEN
               PMAX = DMAX1(PMAX,DABS(DU(I)))
               PULP = DMAX1(PULP,SPACING(U(I)))
            ELSE
               MMAX = DMAX1(MMAX,DABS(DU(I)))
               MULP = DMAX1(MULP,SPACING(U(I)))
            ENDIF
  220    CONTINUE
         WRITE(*,*) 'MAX PRESSURE CORRECTION=',PMAX,
     &              ' MAX PRESSURE ULP=',PULP
         WRITE(*,*) 'MAX FLOW CORRECTION=',MMAX,
     &              ' MAX FLOW ULP=',MULP
         WRITE(*,*) 'LAST TRIAL RESIDUAL=',RNEW
         INFO = 4
         RETURN'''
    original = original.replace(
        '      DOUBLE PRECISION RNORM,RNEW,SNORM,USCALE,LAMBDA',
        '''      DOUBLE PRECISION RNORM,RNEW,SNORM,USCALE,LAMBDA
      DOUBLE PRECISION STERM,FIN,FOUT,QTERM
      DOUBLE PRECISION DENNET,DENTER,ETANET,ETATER
      DOUBLE PRECISION PMAX,MMAX,PULP,MULP
      DOUBLE PRECISION EOS_CZ_RHO''',
        1,
    )
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
