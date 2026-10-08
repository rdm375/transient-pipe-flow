# M8D — Adaptive Integration Accuracy and Cost Validation

Status: CLOSED — benchmark validation milestone
Date: 2026-10-08

## Objective

Evaluate whether adaptive implicit-theta integration reduces
computational cost relative to fixed timestepping at comparable
accuracy in:

1. Pipeline inventory (linepack).
2. Physical inlet/outlet mass-flow imbalance.

Numerical mass-conservation defects are monitored separately.

## Configuration

- Isothermal, single-pipe transient simulator.
- Pipe length: 100 km.
- Pipe diameter: 1 m.
- Spatial intervals: 20 and 40.
- Implicit theta weighting: 0.65.
- Simulation duration: 3600 s.
- Outlet demand: 100 to 105 kg/s over 600 s.
- Inlet pressure: 8 MPa.
- Fixed-step reference: dt = 3.75 s.
- Compiler: gfortran, -O2, -ffixed-line-length-none.
- Timing: median of 11 executable runs.

The same underlying nonlinear transient solver is used by
the fixed and adaptive benchmark drivers.

## Accuracy criteria

Inventory error:

    E_M = max |(M_method(t) - M_ref(t))
               - (M_method(0) - M_ref(0))|

Physical mass-imbalance error:

    E_Q = max |Q_method(t) - Q_ref(t)|

where

    Q(t) = inlet_mass_flow(t) - outlet_mass_flow(t).

The reference trajectory is linearly interpolated at each
method's recorded timestamps.

A fixed-step run qualifies only if both errors are no larger
than the corresponding adaptive-run errors.

The fastest qualifying fixed-step configuration from the
tested sweep is used for the reported speedup.

## Results

| N | Adaptive tolerance | Fixed dt (s) | Speedup |
|---|---:|---:|---:|
| 20 | 0.003 | 6 | 2.407 |
| 20 | 0.010 | 12 | 2.054 |
| 20 | 0.030 | 20 | 1.931 |
| 40 | 0.003 | 6 | 2.372 |
| 40 | 0.010 | 10 | 2.631 |
| 40 | 0.030 | 20 | 2.102 |

Median speedup: approximately 2.24x.

Observed speedup range: 1.93x to 2.63x.

## Findings

1. Adaptive integration reduced measured execution time
   across all six accuracy-matched comparisons.

2. The advantage persisted at both tested spatial
   resolutions.

3. Rejected steps did not eliminate the measured
   performance advantage.

4. Inventory accuracy and physical mass-imbalance accuracy
   must both be considered when comparing integration
   strategies.

5. Small numerical conservation defects do not imply
   small trajectory errors.

## Limitations

1. The 3.75-second reference is a finite-step solution,
   not an exact transient solution.

2. An attempted 1.875-second reference failed at
   t = 2985 s with Newton INFO = 4.
   The failure remains unresolved.

3. Timing includes process startup and CSV output.

4. Newton iteration counts in the adaptive CSV exclude
   iterations performed during rejected attempts.

5. Errors are evaluated at method-specific output times,
   using interpolation of the reference trajectory.

6. The fixed-step comparison uses a discrete timestep
   sweep rather than continuous optimization.

7. Results cover one transient boundary-demand scenario.

8. The adaptive error indicator is heuristic and does
   not certify local or global error.

## Conclusion

Adaptive implicit-theta integration demonstrated a
measured 1.93x–2.63x execution-time advantage over the
fastest accuracy-qualifying fixed-step configurations
in the tested sweep.

The result supports retaining adaptive integration as
a principal execution strategy.

The benchmark does not establish a universal speedup
or a certified global accuracy guarantee.

## Deferred work

- Solver-only performance measurement.
- Total Newton iteration accounting.
- Common observation-grid accuracy comparison.
- Investigation of Newton INFO = 4 near convergence.
- Independent temporal reference convergence.
- Additional boundary-condition scenarios.
- Larger spatial resolutions.

These items are deferred to subsequent validation and
production-hardening work. They are not prerequisites
for closing the present benchmark milestone.

## Reproduction

Run:

    OPENBLAS_NUM_THREADS=1 \
    python3 benchmarks/run_m8d_accuracy_cost.py \
        --n 20 40 \
        --fixed-dt 3.75 4.5 5 6 7.5 9 10 12 15 18 20 24 30 \
        --tol 0.003 0.01 0.03 \
        --duration 3600 \
        --reference-dt 3.75 \
        --repeats 11 \
        --outdir benchmarks/results-m8/m8d-refined

Results:

    benchmarks/results-m8/m8d-refined/summary.csv
