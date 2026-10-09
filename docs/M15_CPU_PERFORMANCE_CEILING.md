# M15 — CPU Performance Ceiling Investigation

**Project:** Isothermal transient single-pipe simulator  
**Status:** PROPOSED — FUTURE WORK  
**Predecessor:** M13 CPU performance parity  
**Independent milestone:** M14 GPU acceleration

## 1. Research objective

Determine the practical CPU performance ceiling for
the validated 100-cell transient pipeline simulation.

The investigation will identify the remaining
computational bottlenecks, quantify their costs,
and evaluate optimizations that exploit the
mathematical structure of the problem.

No performance target in this document is
considered achieved or experimentally established.

## 2. Reference implementations

M13 established the following baselines:

| Implementation | Median CPU time |
|---|---:|
| Fortran77 M13.7 | 90.998 ms |
| C++20 M13.8 | 95.394 ms |

Reference workload:

- 100 spatial cells.
- 15-second timestep.
- 8,000 timesteps.
- 120,000 seconds simulated.
- 1,049 total Newton iterations.

Both implementations use:

    -O3 -march=native -flto

The M13 results are preserved under tag `m13-final`.

## 3. Central observation

The reference simulation performs 1,049 Newton
iterations over 8,000 timesteps.

The average is:

    1049 / 8000 = 0.131125

Newton iterations per timestep.

This suggests that many timesteps converge
without requiring a Newton correction.

However, the aggregate iteration count alone
does not establish the exact distribution
of iterations across timesteps.

The distribution must be measured.

A potentially important remaining cost is
the work required to establish convergence
at the initial Newton guess.

## 4. Performance decomposition

Measure total CPU time attributable to:

1. Initial residual evaluation.
2. Jacobian assembly.
3. Banded linear solution.
4. Newton correction and line search.
5. Physical-property evaluation.
6. State copying and workspace operations.
7. Boundary-condition updates.
8. Linepack and conservation diagnostics.
9. Integration bookkeeping.

Instrumentation should minimize perturbation
of the optimized executable.

Independent kernel benchmarks should be used
to corroborate full-application profiles.

## 5. Newton iteration statistics

Collect the distribution of Newton iterations
per timestep.

Required measurements include:

- Number of zero-iteration timesteps.
- Number of one-iteration timesteps.
- Number requiring multiple iterations.
- Iteration counts during the demand ramp.
- Iteration counts after the ramp.
- Residual norms at initial guesses.
- Residual evaluation counts.
- Jacobian assembly counts.
- Linear solve counts.

The objective is to establish which operations
dominate the common timestep path.

## 6. Optimization hypotheses

### H1 — Specialized zero-correction path

A specialized path may reduce overhead when
the initial Newton guess already satisfies
the convergence criterion.

The convergence test must remain mathematically
equivalent to the validated implementation.

Skipping required residual evaluations is
not an acceptable optimization.

### H2 — Residual specialization

Investigate whether the initial residual can
be simplified further using known relationships
between old-time and initial-guess variables.

Any simplification must preserve the governing
discretization and floating-point behavior
under the chosen validation criteria.

### H3 — Physical-property reuse

Investigate additional invariant quantities
that can be computed once per timestep.

Distinguish carefully between:

- Quantities invariant across Newton iterations.
- Quantities invariant only at the initial guess.
- Quantities dependent on current density.
- Quantities dependent on current mass flow.

### H4 — Specialized nonlinear solution

Investigate predictors, convergence criteria,
and Newton execution paths tailored to the
single-pipe staggered discretization.

A predictor may reduce correction iterations
but increase initial residual cost or invalidate
existing friction-factor reuse opportunities.

Performance must be measured end to end.

### H5 — Specialized linear algebra

Investigate further exploitation of the
fixed Jacobian bandwidth and sparsity pattern.

M13.6 demonstrated that faster isolated
linear algebra does not necessarily yield
substantial full-simulation improvements.

The potential gain must be bounded by the
fraction of execution time spent in linear solves.

### H6 — Reduced timestep overhead

Investigate:

- State copying.
- Array traversal.
- Repeated parameter validation.
- Boundary-condition evaluation.
- History storage.
- Conservation diagnostics.

Preserve the externally observable results
and required diagnostics.

## 7. Performance modeling

Use measured runtime fractions to construct
Amdahl-law bounds.

For an operation consuming fraction f of
total execution time, a speedup s produces:

    S = 1 / ((1 - f) + f/s)

An optimization's theoretical upper bound
must be calculated before substantial
implementation effort.

Investigate whether multiple improvements
interact or eliminate overlapping work.

## 8. Exploratory performance targets

The following ranges are hypotheses,
not validated predictions.

| Target | Total CPU time | Status |
|---|---:|---|
| M13 baseline | ~91 ms | Measured |
| Incremental improvement | 60–80 ms | Exploratory |
| Aggressive specialization | 20–40 ms | Speculative |
| Alternative algorithm | Below 20 ms | Highly speculative |

The profiling results may invalidate
any or all of these targets.

The objective is to establish an
evidence-based practical ceiling,
not to force a predetermined speedup.

## 9. Numerical acceptance criteria

All optimizations must preserve the
validated physical model and discretization
unless explicitly introduced as a
separate numerical-method experiment.

Validation should include:

- Residual parity.
- Jacobian parity.
- Newton solver parity.
- Full transient integration parity.
- Mass-conservation diagnostics.
- Difficult transient cases.
- Multiple spatial resolutions.
- Multiple timestep sizes.

Numerical differences must be investigated
causally rather than concealed by
arbitrarily widening tolerances.

Changes to convergence criteria or
floating-point operation ordering require
explicit documentation and validation.

## 10. Experimental methodology

Use the M13 benchmark methodology as a baseline:

- Matched compiler optimization settings.
- Preserved reference executables.
- External process CPU-time measurements.
- Warm-up runs.
- Alternating execution order.
- Multiple repetitions.
- Median and distribution reporting.
- Executable checksums.
- Numerical output comparison.

Where practical, add confidence intervals,
CPU-frequency observations, and thermal-state
measurements.

Microbenchmark improvements must be verified
against total simulation runtime.

## 11. Proposed execution phases

### M15.1 — Profiling and cost decomposition

Produce a measured breakdown of the
remaining CPU execution time.

### M15.2 — Newton iteration analysis

Characterize zero-correction timesteps
and the cost of initial convergence checks.

### M15.3 — Residual and property specialization

Implement and validate the most promising
mathematically justified optimizations.

### M15.4 — Linear algebra and timestep overhead

Evaluate remaining opportunities using
measured Amdahl-law bounds.

### M15.5 — Final performance ceiling assessment

Compare the optimized implementations
against the M13 reference baseline.

Document:

- Achieved speedups.
- Remaining dominant costs.
- Numerical validation.
- Practical lower-bound estimates.
- Workload dependence.
- Recommendations for future research.

## 12. Relationship to M14

M14 investigates GPU acceleration.

M15 investigates the CPU performance ceiling.

These are independent research directions
sharing the validated M13 numerical baseline.

Findings from M15 may inform GPU kernel design,
particularly regarding residual evaluation,
nonlinear iteration structure, and redundant
physical-property calculations.

**M15 STATUS: PROPOSED.**
