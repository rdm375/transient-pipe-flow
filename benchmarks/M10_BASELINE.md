# M10.1 CPU performance baseline

The benchmark calls the existing M9E `INTEGRATE_TRANSIENT` routine.

Run:

    make m10-baseline

Or:

    python3 benchmarks/run_m10_baseline.py --repetitions 7

The sweep varies spatial resolution and timestep size while holding
the physical pipe length and simulation duration fixed.

Outputs:

- `results-m10/baseline.csv`: individual observations.
- `results-m10/baseline.json`: provenance.

`cpu_seconds` measures Fortran `CPU_TIME` around the integrator.
`wall_seconds` includes process startup and output handling.

Numerical diagnostics include final outlet pressure, final linepack,
total Newton iterations, and maximum mass-conservation defect.

The mass-conservation defect is numerical error, not the physical
inlet-minus-outlet flow imbalance.

This is an end-to-end baseline. Component-level profiling is deferred
to M10.2. The current solver supports at most 100 spatial intervals.
