# M10 CPU Performance and Numerical Parity

## Methodology

- Solver: isothermal transient single-pipe simulation.
- Dense reference: original dense Newton linear solver.
- Optimized backend: direct-banded Jacobian and linear solver.
- Compiler: GNU Fortran 15.2.0.
- Optimization: `-O3 -march=native`.
- Simulation duration: 1,200 seconds.
- Mesh intervals: 10, 20, 40, 80, 100.
- Timestep sizes: 15, 30, 60 seconds.
- Repetitions: five per configuration and backend.
- Timing: Fortran `CPU_TIME` around `INTEGRATE_TRANSIENT`.
- Reported execution times: median of five repetitions.
- Historical run revision: `9e1c57d` (uncommitted M10.4 development state).

## Performance

| N | dt (s) | Dense (ms) | Banded (ms) | Speedup |
|---:|---:|---:|---:|---:|
| 10 | 15 | 1.3690 | 0.9210 | 1.49x |
| 10 | 30 | 0.7360 | 0.4630 | 1.59x |
| 10 | 60 | 0.3800 | 0.2420 | 1.57x |
| 20 | 15 | 4.6520 | 1.7330 | 2.68x |
| 20 | 30 | 2.3860 | 0.8720 | 2.74x |
| 20 | 60 | 1.3090 | 0.4440 | 2.95x |
| 40 | 15 | 23.0170 | 3.2750 | 7.03x |
| 40 | 30 | 11.7040 | 1.6540 | 7.08x |
| 40 | 60 | 5.9140 | 0.8400 | 7.04x |
| 80 | 15 | 162.5130 | 6.5470 | 24.82x |
| 80 | 30 | 76.5040 | 3.2860 | 23.28x |
| 80 | 60 | 42.1350 | 1.6230 | 25.96x |
| 100 | 15 | 348.1550 | 8.8820 | 39.20x |
| 100 | 30 | 170.8990 | 4.1010 | 41.67x |
| 100 | 60 | 109.0230 | 2.0670 | 52.74x |

## Numerical parity

- Configurations: 15.
- Individual comparisons: 75.
- Newton iteration mismatches: 0.
- Maximum `outlet_pressure` difference: 0.000000000000e+00.
- Maximum `linepack` difference: 0.000000000000e+00.
- Maximum `max_mass_defect` difference: 0.000000000000e+00.

The differences above compare the values recorded in the benchmark CSV files. Zero differences demonstrate agreement of the recorded results, not bitwise equality of every intermediate floating-point operation.

## Provenance

The original benchmark metadata is preserved in:

- `baseline.json` (M10.1 reference).
- `comparison-dense.json`.
- `comparison-banded.json`.

The comparison JSON files record the original working-tree state, compiler, platform, benchmark configuration, and executable paths.
