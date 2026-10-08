# Isothermal Pipe Simulator — Development Roadmap

**Status:** Proposed development roadmap
**Updated:** 2026-10-08
**Current branch:** `development/m9-production-hardening`

## 1. Project Direction

The project is developing a scientifically defensible, open-source
transient pipeline simulator.

The immediate objective is a reliable, validated Fortran reference
implementation for a single isothermal pipe.

Subsequent development may include:

- Reproducible CPU performance characterization.
- A numerically equivalent modern C++ implementation.
- CPU optimization using the structure of the discrete equations.
- GPU feasibility studies and implementation where justified.
- Cross-platform numerical validation.
- Adaptive GPU execution.
- Extended physical models.
- A reproducible research-grade release.

The single-pipe solver should remain independent of pipeline-network
topology.

Future network simulation should build upon validated pipe components
rather than complicating the present single-pipe formulation.

## 2. Milestone Overview

| Milestone | Focus | Primary deliverable |
|---|---|---|
| M9 | Production hardening | Reliable, validated Fortran solver |
| M10 | Performance characterization | Reproducible CPU baseline |
| M11 | C++ implementation | Numerically equivalent C++ solver |
| M12 | CPU optimization | Optimized solver execution |
| M13 | GPU architecture | Feasibility benchmarks and design |
| M14 | GPU implementation | Working GPU transient solver |
| M15 | CPU/GPU validation | Explained numerical equivalence |
| M16 | Adaptive GPU execution | Efficient adaptive GPU simulations |
| M17 | Extended physical models | Independently validated extensions |
| M18 | Research-grade release | Reproducible open-source release |

These milestones describe a proposed direction rather than fixed
implementation commitments.

GPU development is contingent upon measured feasibility.

## 3. M9 — Production Hardening

**Objective:** Make the existing Fortran implementation reliable under
valid, invalid, and numerically difficult operating conditions.

### M9A — Newton Solver Reliability

Completed:

- Define iteration counts on Newton failure paths.
- Investigate the observed line-search stagnation.
- Demonstrate guarded convergence experimentally.
- Document potential floating-point refinements.

Remaining:

- Verify deterministic status reporting for every Newton exit path.
- Test singular and poorly conditioned linear systems.
- Test line-search failure.
- Test maximum-iteration exhaustion.
- Ensure unsuccessful solves do not corrupt accepted simulation states.

The guarded stagnation rule remains experimental.

Potential arithmetic reformulations are documented separately in
`docs/future-numerical-refinements.md`.

**Exit criterion:** Every Newton outcome is explicit, deterministic,
and covered by regression tests.

### M9B — Physical Input Validation

Validate:

- Pipe length and diameter.
- Temperature.
- Timestep and simulation duration.
- Gas constants and compressibility factors.
- Pressure and density.
- Friction-model parameters.
- Boundary-condition consistency.
- Finite numerical inputs.

Test:

- Zero flow.
- Reverse flow.
- Near-zero flow.
- Invalid physical states.
- Nonfinite values.
- Degenerate configurations.

**Exit criterion:** Invalid inputs produce predictable failure
statuses without undefined behavior or uncontrolled floating-point
exceptions.

### M9C — Integration Failure Handling

Establish transactional timestep semantics.

A timestep must either:

1. Complete successfully and commit its new state.
2. Fail without changing the previously accepted state.

Distinguish:

- Recoverable adaptive timestep rejection.
- Nonlinear solver failure.
- Invalid physical state.
- Unrecoverable integration failure.

Test minimum and maximum timestep behavior.

**Exit criterion:** Failed and rejected timesteps cannot corrupt the
accepted trajectory.

### M9D — Numerical Regression Coverage

Extend the existing regression suite to cover:

- Steady-state preservation.
- Forward and reverse flow.
- Zero-flow equilibrium.
- Pressure transients.
- Demand transients.
- Spatial refinement.
- Temporal refinement.
- Global mass conservation.
- Analytic Jacobian consistency.
- Compiler optimization sensitivity.
- Deterministic failure behavior.

Preserve justified numerical tolerances.

Investigate discrepancies rather than widening tolerances without
explanation.

**Exit criterion:** The documented regression suite covers the
supported physical and numerical operating domain.

### M9E — Reproducible Solver Interface

Define stable interfaces for:

- Physical parameters.
- Initial conditions.
- Boundary schedules.
- Solver configuration.
- Integration results.
- Diagnostics.
- Failure statuses.

Separate the numerical engine from test and benchmark drivers.

Avoid unnecessary abstractions in the Fortran 77 reference
implementation.

**Exit criterion:** Independent programs can invoke the solver
consistently.

### M9F — Production Release Candidate

Complete:

- Clean debug and optimized builds.
- Compiler-warning review.
- Regression testing.
- Reproducible reference simulations.
- Documentation of supported operating conditions.
- Documentation of known limitations.
- Versioned release candidate.

**Exit criterion:** A tagged, reproducible Fortran reference solver.

## 4. M10 — CPU Performance Baseline

**Objective:** Establish trustworthy performance measurements before
optimization or porting.

Measure performance versus:

- Spatial resolution.
- Timestep.
- Requested accuracy.
- Newton iteration count.
- Boundary-condition complexity.

Separate costs for:

- Residual evaluation.
- Jacobian assembly.
- Linear solution.
- Newton iteration.
- Timestep control.

Use the M8D adaptive benchmark as an initial reference.

Its approximately 2.24x adaptive speedup applies to the qualifying
benchmark configuration, not necessarily other workloads.

**Exit criterion:** Reproducible CPU benchmarks relating runtime
to numerical accuracy.

## 5. M11 — Modern C++ Implementation

**Objective:** Develop a C++ implementation numerically equivalent
to the validated Fortran reference.

This need not be a line-by-line translation.

Preserve the mathematical discretization while introducing clear
interfaces for states, physical models, boundary conditions, and
integration.

### M11A — EOS and Friction Parity

Compare EOS and friction calculations.

### M11B — Steady Residual Parity

Compare steady-state residual evaluations.

### M11C — Transient Residual Parity

Compare transient residual evaluations.

### M11D — Jacobian Parity

Compare analytic Jacobian entries and sparsity.

### M11E — Newton Parity

Compare corrections, iteration behavior, and failure statuses.

### M11F — Timestep Parity

Compare individual implicit timesteps.

### M11G — Trajectory Parity

Compare complete transient trajectories.

### M11H — Adaptive Integration Parity

Compare adaptive timestep decisions and conservation diagnostics.

**Exit criterion:** The C++ implementation agrees with the Fortran
reference within justified numerical tolerances.

## 6. M12 — CPU Optimization

**Objective:** Improve CPU performance without sacrificing
numerical correctness.

Investigate the structure of the single-pipe Jacobian.

Potential approaches:

- Banded direct solution.
- Block elimination.
- Sparse factorization.
- Reuse of matrix structure.
- Reduced assembly overhead.

Establish bandwidth, operation counts, and profiling evidence before
selecting an algorithm.

**Exit criterion:** Measured performance improvements at equivalent
physical accuracy.

## 7. M13 — GPU Architecture and Feasibility

**Objective:** Determine which workloads benefit from GPU execution.

Evaluate:

1. A single pipe simulated on one GPU.
2. Many independent pipes simulated concurrently.
3. Batched parameter sweeps and uncertainty studies.

Investigate the implications of implicit timestep dependencies and
linear-system solution.

Compare against the optimized CPU implementation.

**Exit criterion:** A measured GPU execution strategy with a credible
performance advantage for a defined workload.

A single small pipe may not justify GPU execution.

## 8. M14 — GPU Implementation

**Objective:** Implement the GPU strategy selected in M13.

Potential components:

- GPU residual evaluation.
- GPU Jacobian assembly.
- Structured linear-system solution.
- Batched Newton iteration.
- Device-resident state management.

Minimize unnecessary host/device transfers.

**Exit criterion:** Complete transient simulations executing on
the GPU.

## 9. M15 — CPU/GPU Numerical Validation

**Objective:** Establish numerical agreement between the CPU
reference and GPU implementations.

Compare:

- Residuals.
- Jacobians.
- Newton corrections and convergence.
- Pressure trajectories.
- Mass-flow trajectories.
- Inventory.
- Physical mass imbalance.
- Numerical conservation defects.
- Adaptive timestep decisions.

Use the most precise available physical and numerical constants.

Investigate material discrepancies rather than widening tolerances
to conceal them.

**Exit criterion:** Material CPU/GPU numerical differences are
explained or eliminated.

## 10. M16 — Adaptive GPU Execution

**Objective:** Execute adaptive transient simulations efficiently
on GPU hardware.

Investigate:

- Timestep bucketing.
- Grouping by convergence behavior.
- GPU-resident timestep control.
- Rejection and retry.
- Workload scheduling.
- Divergent simulation trajectories.

**Exit criterion:** Measured advantage for relevant adaptive
GPU workloads.

## 11. M17 — Extended Physical Models

**Objective:** Extend physical capabilities after the numerical
infrastructure is stable.

Potential additions:

- Temperature-dependent gas properties.
- More sophisticated equations of state.
- Alternative friction correlations.
- Additional boundary-condition models.
- Non-isothermal transient flow.

Non-isothermal flow requires a substantial mathematical extension
and must not be treated as merely replacing the EOS.

Each extension requires independent validation.

**Exit criterion:** New physical models have documented validation
cases and conservation tests.

## 12. M18 — Research-Grade Release

**Objective:** Consolidate the validated implementations into a
reproducible open-source scientific-computing project.

Potential deliverables:

- Fortran reference implementation.
- Validated C++ implementation.
- GPU implementation where justified.
- Benchmark datasets.
- Numerical-method documentation.
- Validation reports.
- Reproducible builds.
- Clearly stated limitations.
- Technical report or research publication.

**Exit criterion:** A documented, reproducible release suitable
for independent scientific evaluation.

## 13. Dependencies

The principal development path is:

    M9  Production hardening
     |
     +-- M10 CPU performance baseline
     |    |
     |    +-- M12 CPU optimization
     |         |
     |         +-- M13 GPU feasibility
     |              |
     |              +-- M14 GPU implementation
     |                   |
     |                   +-- M15 CPU/GPU validation
     |                        |
     |                        +-- M16 Adaptive GPU execution
     |
     +-- M11 C++ implementation
     |    |
     |    +-- M12 CPU optimization
     |
     +-- M17 Extended physics
     |
     +-- M18 Research-grade release

M10 and M11 may proceed partly in parallel.

M12 should use the validated C++ implementation and the M10
performance baseline.

M13 is a feasibility gate, not an unconditional commitment
to GPU development.

M17 may proceed independently once its required numerical
infrastructure is sufficiently mature.

M18 consolidates completed, validated capabilities.

## 14. Immediate Development Sequence

The recommended sequence is:

    M9A -> M9B -> M9C -> M9D -> M9E -> M9F

followed by:

    M10 -> M11 -> M12 -> M13

The current priority is completing M9 production hardening.

The M9 floating-point cancellation investigation is deferred.
Its findings and experimental diagnostics are preserved in the
repository.

## 15. Engineering Principles

Throughout development:

1. Preserve the Fortran implementation as a numerical reference.
2. Keep the single-pipe model independent of network topology.
3. Validate mathematical correctness before optimization.
4. Maintain global mass-conservation diagnostics.
5. Distinguish physical imbalance from numerical conservation error.
6. Preserve reproducible benchmarks.
7. Investigate numerical discrepancies causally.
8. Avoid unjustified relaxation of tolerances.
9. Require evidence before adopting new algorithms or hardware.
10. Keep experimental features separate from production defaults.
