# Adaptive transient integration

## M8C implementation

The adaptive integrator advances the staggered, isothermal
single-pipe model using the implicit theta method.

The implementation provides:

- Adaptive timesteps without step doubling.
- A history-based temporal error indicator.
- Separate absolute and relative error scales for pressure,
  mass flow, linepack, and boundary mass imbalance.
- Continuous piecewise-linear inlet-pressure schedules.
- Continuous piecewise-linear outlet-mass-flow schedules.
- Exact timestep termination at boundary-schedule knots.
- Optional acoustic CFL accuracy restriction.
- Rejected-step retries.
- Linepack and mass-conservation diagnostics.

The existing fixed-step integrator remains available.

## Error estimation

The timestep controller uses differences between successive
discrete solution slopes to estimate temporal variation.

This indicator is a heuristic approximation to local
temporal truncation error.

It is not a rigorous upper bound on local or global error.

The controller tolerance therefore must not be interpreted
as a guaranteed bound on pressure, mass flow, or linepack
trajectory error.

Accuracy must be established independently using
convergence studies and reference calculations.

## Mass conservation

The staggered discretization satisfies the discrete
inventory balance

    M[n+1] - M[n]
      = dt * ((1-theta) * (min[n]-mout[n])
              + theta * (min[n+1]-mout[n+1]))

up to nonlinear-solver and floating-point errors.

Two different quantities must be distinguished:

1. Physical mass imbalance:
   inlet mass flow minus outlet mass flow.

2. Numerical conservation defect:
   the discrepancy between inventory change and
   integrated boundary mass imbalance.

A small numerical conservation defect does not imply
a small inventory trajectory error.

## Boundary schedules

M8C supports continuous piecewise-linear schedules.

Schedule times must be strictly increasing, start at zero,
and cover the requested simulation duration.

Timesteps terminate at schedule knots.

Discontinuous boundary changes are not yet supported by
the schedule representation.

## Current limitations

- Maximum supported spatial resolution is 100 intervals.
- Boundary schedules are continuous and piecewise linear.
- Error estimation is heuristic rather than certified.
- Global inventory accuracy is not guaranteed by ETOL.
- Newton work accounting depends on solver return values.
- Production failure-path coverage remains incomplete.

## Validation

The M8C regression suite includes:

- A variable outlet-demand schedule.
- A variable inlet-pressure and outlet-demand schedule.
- Boundary-knot alignment.
- Final-time completion.
- Accepted-step indicator checks.
- Integrated inventory conservation.
- Attempt-count consistency.

These tests establish regression behavior, not general
accuracy guarantees.

## Next validation milestone

Compare adaptive integration against fixed-step integration
at matched maximum inventory trajectory error and matched
maximum physical mass-imbalance error.

Record:

- Accepted timesteps.
- Rejected timesteps.
- Total Newton iterations.
- Wall-clock execution time.
- Maximum inventory trajectory error.
- Maximum mass-imbalance trajectory error.
- Maximum numerical conservation defect.

Use an independently converged temporal reference.
