# M7b: wave-response instrumentation (first increment)

Run `make wave-analysis` to generate `benchmarks/results-m7b/sensors.csv`, `metrics.csv`, and `REPORT.md`. This does not modify the production solver. The harness uses the same constant-Z model and demand histories as M7a.

54 combinations: forcing=1 (step), 2 (pulse), 3 (smooth cosine ramp); N=20,40; dt=15,60,120 seconds; theta=0.50,0.65,1.00; duration 1800 seconds. Pressure sensors are at nodal indices N/4, N/2, 3N/4, N. All quantities are in SI units.

`step_mass_defect_kg = (linepack_new - linepack_old) - dt * [theta * (min_new - mout_new) + (1-theta) * (min_old - mout_old)]`. `cum_mass_defect_kg` compares total linepack change with accumulated theta-weighted boundary mass flow.

`metrics.csv` contains sampled extrema, their times, a 10%-of-local-maximum threshold crossing, and slope-reversal counts. Arrival times are sampling-dependent and cannot by themselves establish characteristic propagation speed. A step or pulse sampled at endpoints has boundary forcing quadrature effects. Slope reversals are a descriptive heuristic, not an oscillation stability proof. The linearized dispersion relation and independently validated damping predictions are still open tasks.

Run `make check` to verify the existing M1–M5 suite. The new diagnostics assert cumulative conservation defect <1e-7 relative to initial linepack for the enumerated cases. This is a regression threshold, not an error bound.

## Linearized modal analysis (second increment)

Run `make modal-analysis` (requires NumPy) to generate `modal_matrix.csv`,
`modal_metrics.csv`, and `MODAL_REPORT.md` in `results-m7b/`. The F77
matrix exporter calls the *production analytical Jacobian* at the exact
steady state and extracts the differential subsystem after imposing fixed
inlet pressure and eliminating algebraic inlet flow. Each eigenvalue is
mapped through the exact theta-method rational amplification function.
This predicts discrete linearized decay and sampled frequency, but does
not itself validate the nonlinear sensor traces or a continuum wave speed.

## M7b small-perturbation closure experiment

Run `make small-verify` after installing NumPy. This target runs the
54-case sensor suite, extracts the production-Jacobian modal operator,
and compares nonlinear Fortran time histories against a linearized
forced theta-method simulation for a smooth outlet-demand perturbation.

The experiment uses amplitudes 1, 0.5, 0.25, 0.125 kg/s; grids N=40;
timesteps 15 and 60 s; and theta weights 0.5, 0.65 and 1.0. The
linearized forcing derivative is the negative inverse of the outlet
half-cell mass coefficient. It retains the production inlet-pressure
constraint and the full staggered discrete operator.

Outputs: `benchmarks/results-m7b/small_nonlinear.csv`,
`small_metrics.csv`, and `SMALL_REPORT.md`. The verification criterion
is observed near-quadratic pressure discrepancy under amplitude
refinement (minimum observed order 1.5 over the tested factor-of-eight
range). This tests the consistency of the linearized, forced,
time-discrete model, **not** physical wave speed, global temporal error,
or individual eigenmode excitation. The modal amplification formula
continues to provide an exact linear-theta prediction for each
semidiscrete eigenvalue.
