#!/usr/bin/env python3
"""Compile and run M8 adaptive Fortran driver; compare against M7c-18A fixed-step reference."""
import argparse
import csv
import subprocess
from pathlib import Path
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
def read(path):
 with path.open(newline='') as f: rows=list(csv.DictReader(f))
 return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]}

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--n',type=int,nargs='+',default=[20,40])
 p.add_argument('--tol',type=float,nargs='+',default=[0.003,0.01,0.03])
 p.add_argument('--duration',type=float,default=3600)
 p.add_argument('--hmax',type=float,default=120)
 p.add_argument('--reference-dt',type=float,default=3.75)
 p.add_argument('--outdir',type=Path,default=Path('benchmarks/results-m8/adaptive-inventory'))
 a=p.parse_args()
 src=Path('benchmarks/adaptive_inventory_m8.f')
 for f in [src,*map(Path,SOURCES),Path('benchmarks/inventory_history_m7c18a.f')]:
  if not f.is_file(): p.error(f'Missing file: {f}')
 a.outdir.mkdir(parents=True,exist_ok=True)
 exe=a.outdir/'adaptive_inventory_m8'
 fixed=a.outdir/'inventory_history_m7c18a'
 flags=['gfortran','-O2','-ffixed-line-length-none']
 subprocess.run([*flags,'-o',str(exe),*SOURCES,str(src)],check=True)
 subprocess.run([*flags,'-o',str(fixed),*SOURCES,'benchmarks/inventory_history_m7c18a.f'],check=True)
 summary=[]
 for n in a.n:
  refpath=a.outdir/f'n{n}_reference_dt{a.reference_dt:g}.csv'
  subprocess.run([str(fixed),str(a.reference_dt),str(a.duration),str(refpath),str(n)],check=True)
  ref=read(refpath)
  for tol in a.tol:
   path=a.outdir/f'n{n}_tol{tol:g}.csv'
   subprocess.run([str(exe),str(path),str(n),str(a.duration),str(tol),str(a.hmax)],check=True)
   r=read(path)
   t=r['time_s']; delta=r['inventory_kg']-np.interp(t,ref['time_s'],ref['inventory_kg'])
   dq=r['imbalance_kg_s']-np.interp(t,ref['time_s'],ref['imbalance_kg_s'])
   result=dict(n=n,tol=tol,steps=len(t)-1,rejections=int(np.sum(r['rejects_before_accept'])),
     newton_iterations=int(np.sum(r['newton_iterations'])),
     min_dt_s=float(np.min(r['dt_s'][1:])),max_dt_s=float(np.max(r['dt_s'][1:])),
     max_inventory_change_error_kg=float(np.max(np.abs(delta-delta[0]))),
     final_inventory_change_error_kg=float(delta[-1]-delta[0]),
     max_imbalance_error_kg_s=float(np.max(np.abs(dq))),
     max_conservation_defect_kg=float(np.max(np.abs(r['step_balance_defect_kg'][1:]))))
   summary.append(result)
   print(f'n={n} tol={tol:g} steps={result["steps"]} rejected={result["rejections"]} '
         f'max inventory error={result["max_inventory_change_error_kg"]:.5g} kg '
         f'max imbalance error={result["max_imbalance_error_kg_s"]:.5g} kg/s',flush=True)
 with (a.outdir/'summary.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
 print('Wrote',a.outdir/'summary.csv')
 print('NOTE: reference is a finite-step solution; error estimator is heuristic, not certified.')
if __name__=='__main__': main()
