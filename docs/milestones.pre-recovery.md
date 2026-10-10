# Isothermal Pipe Simulator Milestones

## Current development roadmap — CPU to GPU

This roadmap supersedes the milestone numbering in the historical
planning sections below. Earlier research objectives remain relevant,
but are not automatically considered completed or assigned to the
new milestones.

The development sequence is:

1. Optimize and freeze the Fortran77 CPU reference.
2. Port the validated numerical algorithm to C++.
3. Establish C++ numerical and performance parity.
4. Introduce GPU acceleration.
5. Optimize and validate CPU/GPU execution.

All implementations retain the same governing equations, conservative
spatial discretization, theta-method integration, boundary conditions,
and numerical-assurance requirements unless a change is explicitly
specified and independently validated.

### M10 — Banded CPU solver optimization

**Status: COMPLETE — published on `development/m10-performance`.**

- [x] M10.1: Establish optimized dense CPU baseline.
- [x] M10.2: Profile the dense implementation.
- [x] M10.3: Establish Jacobian bandwidth and solver parity.
- [x] M10.4: Implement direct-banded Jacobian assembly and Newton solve.
- [x] M10.5: Select banded Newton as the default optimized CPU backend.
- [x] Preserve the dense implementation for independent comparison.
- [x] Verify dense/banded numerical parity across 75 benchmark runs.
- [x] Preserve benchmark measurements and provenance.
- [x] Publish implementation and performance evidence.

Reference workload: N=100, dt=15 s, duration=1200 s,
80 timesteps and 160 Newton iterations.

Measured dense/banded CPU times: approximately 330 ms and 8.3 ms
in the final direct comparison, respectively.

The full benchmark sweep is documented in
`benchmarks/results-m10/M10_PERFORMANCE.md`.

### M11 — Fortran77 performance optimization

**Status: PLANNED — next milestone.**

- [ ] M11.1: Reprofile the optimized banded solver.
- [ ] M11.2: Investigate residual and Jacobian evaluation costs.
- [ ] M11.3: Investigate banded factorization and substitution costs.
- [ ] M11.4: Investigate memory access, redundant work, and compiler
      optimization behavior.
- [ ] M11.5: Freeze the optimized F77 reference and benchmark suite.

Requirements:

- Profile before optimizing.
- Measure changes against the frozen M10 baseline.
- Preserve Newton convergence and failure contracts.
- Preserve physical and numerical conservation diagnostics.
- Test multiple mesh sizes, timesteps, and flow directions.
- Reject unexplained numerical discrepancies.
- Record compiler flags, hardware, revision, and benchmark variability.
- Avoid relaxed floating-point semantics in the reference configuration.

### M12 — C++20 numerical port

**Status: PLANNED.**

- [ ] Implement the same mathematical model in C++20.
- [ ] Preserve unknown ordering and banded Jacobian structure initially.
- [ ] Implement equivalent residual, Jacobian, Newton, and time integration.
- [ ] Preserve boundary half-cell conservation equations.
- [ ] Match physical outputs and numerical diagnostics against F77.
- [ ] Verify convergence, failure behavior, and reversed-flow cases.
- [ ] Establish a cross-language regression suite.
- [ ] Freeze the numerically validated C++ baseline.

Initially prioritize numerical transparency over architectural redesign.
Keep the F77 reference permanently available.

### M13 — C++ CPU performance parity

**Status: PLANNED.**

- [ ] Establish comparable F77 and C++ benchmark builds.
- [ ] Profile C++ before optimizing.
- [ ] Investigate data layout, aliasing, inlining, and vectorization.
- [ ] Compare GCC and Clang where available.
- [ ] Establish C++ performance within 5% of optimized F77,
      or faster, on representative workloads.
- [ ] Assess measurement variability before accepting small differences.
- [ ] Preserve cross-language numerical and conservation parity.
- [ ] Freeze both CPU implementations and their benchmark evidence.

The 5% performance target is an engineering objective, not a
numerical-correctness tolerance.

### M14 — CUDA implementation

**Status: PLANNED.**

- [ ] Establish the GPU execution and memory architecture.
- [ ] Implement FP64 numerical kernels initially.
- [ ] Validate GPU residual and Jacobian evaluation against CPU.
- [ ] Investigate structured GPU linear solvers.
- [ ] Implement end-to-end transient integration.
- [ ] Investigate single-pipe execution.
- [ ] Investigate batched independent-pipe execution.
- [ ] Preserve physical models, discretization, and failure semantics.

Single-pipe latency and batched throughput must be measured separately.

### M15 — GPU performance optimization and validation

**Status: PLANNED.**

- [ ] Establish CPU/GPU numerical parity over a validation matrix.
- [ ] Verify Newton convergence and failure handling.
- [ ] Verify linepack and global mass conservation.
- [ ] Distinguish physical inlet/outlet mass imbalance from numerical
      conservation defect.
- [ ] Profile GPU kernels and data transfers.
- [ ] Optimize memory layout, occupancy, and solver execution.
- [ ] Measure single-simulation latency and batched throughput.
- [ ] Determine CPU/GPU crossover regimes.
- [ ] Preserve reproducible performance and validation evidence.

### Deferred research objectives

The following objectives from the earlier roadmap remain open and
are not implicitly completed by the CPU-to-GPU milestones:

- Versioned numerical-assurance evidence and explainable reports.
- Interactive numerical exploration and visualization.
- Modular EOS, friction, momentum, and boundary-condition models.
- Peng-Robinson EOS integration and extended physical-model validation.
- Adaptive surrogate modeling and TDAR parameter-space research.
- Adjoint sensitivities, state estimation, and network extensions.

These should receive their own milestones when their scope is adopted.

---

## Historical specification and milestone plan

The sections below preserve the original specification, validation
history, and research roadmap. Some milestone numbers below refer
to an earlier planning sequence and are superseded by the current
CPU-to-GPU roadmap above.


## V0 — Mathematical and numerical specification

- [x] Define isothermal single-pipe governing equations.
- [x] Define reduced and full momentum models.
- [x] Define replaceable EOS and friction closures.
- [x] Select staggered conservative spatial discretization.
- [x] Define theta-method temporal discretization.
- [x] Derive fully discrete interior residuals.
- [x] Derive conservative boundary half-control-volume equations.
- [x] Establish exact discrete global mass balance.
- [x] Derive analytic constant-Z, constant-friction steady solution.
- [x] Establish exact nodal steady-state p-squared property.
- [x] Define validation ladder and public-source provenance.

## V1 — F77 reference implementation

### M1 — Steady residual oracle

- [x] Implement constant-Z analytic steady solution.
- [x] Implement constant-friction reduced momentum residual.
- [x] Implement boundary continuity residuals.
- [x] Verify exact-zero steady continuity residual.
- [x] Verify roundoff-level steady momentum residual.
- [x] Sweep spatial resolutions from N=1 through N=1000.
- [x] Add scale-aware relative momentum residual criterion.

### M2 — Residual and Jacobian kernels

- [x] Refactor residual evaluation into reusable F77-style routines.
- [x] Define the complete transient unknown-vector ordering.
- [x] Implement complete theta-method transient residual.
- [x] Derive analytic Jacobian entries.
- [x] Implement analytic Jacobian assembly.
- [x] Implement independent finite-difference Jacobian.
- [x] Compare analytic and finite-difference Jacobians.
- [x] Test Jacobian at non-equilibrium states.
- [x] Verify expected Jacobian sparsity/locality.

### M3 — Newton nonlinear solve

- [x] Implement linear solve for the reference Jacobian.
- [x] Implement Newton iteration.
- [x] Implement residual-based convergence criteria.
- [x] Implement damped Newton/backtracking.
- [x] Verify recovery of analytic steady state from perturbed guesses.
- [x] Verify failure reporting for nonconvergent cases.

### M4 — Single time step

- [x] Implement one theta-method transient step.
- [x] Preserve an exact steady state for one step.
- [x] Verify per-step discrete global mass balance.
- [x] Verify boundary histories at t^n and t^(n+1).
- [x] Verify inlet boundary flow diagnostic.

### M5 — Transient integration

- [x] Implement multi-step integration.
- [x] Implement downstream mass-flow ramp case.
- [x] Record pressure, mass-flow, linepack, and boundary histories.
- [x] Verify global mass conservation over complete simulation.
- [x] Verify relaxation to the analytic steady solution.

### M6 — Numerical characterization

- [x] Temporal refinement study.
- [x] Spatial refinement study.
- [x] Theta study: 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 1.00.
- [x] Characterize numerical damping and transient accuracy.
- [x] Establish reference cases and tolerances.

### M7a — Transient disturbance characterization

- [x] Implement downstream demand step, rectangular pulse, and smooth ramp cases.
- [x] Compare seven theta values against a fine-timestep reference.
- [x] Export pressure, inlet-flow, and linepack histories as CSV.
- [x] Exercise M7a characterization and M1–M6 checks on Dell 7710.
- [x] Freeze and commit M7a (`192fe21`).

### M7b — Numerical damping and wave-response characterization

- [x] Export sampled pressure extrema and timing; flow extrema remain.
- [x] Add preliminary pressure slope-reversal heuristic; oscillation decay and classification remain.
- [x] Derive exact theta amplification of the linearized constrained semidiscrete operator from the production Jacobian; continuum wave-speed derivation remains.
- [x] Tabulate predicted discrete modal damping and aliased frequency alongside the nonlinear sensor evidence; direct modal-excitation agreement remains to be tested.
- [x] Export four spatial pressure sensors and threshold-crossing arrivals; phase and attenuation interpretation remain.
- [x] Generate multi-resolution sensor histories (N=20,40; dt=15,60,120 s); rigorous reference-error separation remains.
- [x] Exercise dt=120 s and check linear modal stability; nonlinear solver failure envelopes and larger-step limits remain open.
- [x] Record limitations and diagnostic definitions in generated M7b report.

### M8 — Numerical assurance framework

- [ ] Define versioned, solver-independent evidence schema and units.
- [ ] Export per-step conservation residual and cumulative mass-balance defect.
- [ ] Export nonlinear residual histories, iteration counts, line-search factors,
      convergence reasons, and failure states.
- [ ] Define independent time/space refinement studies and error estimates.
- [ ] Check positivity and EOS/friction correlation applicability domains.
- [ ] Attach configuration, source revision, numerical methods, compiler,
      tolerances, and reference provenance to each run.
- [ ] Classify each finding as proven, verified, estimated, or unassessed.
- [ ] Explicitly distinguish numerical assurance from physical validation.
- [ ] Fail closed on missing evidence; never collapse warnings into PASS.
- [ ] Validate evidence records and test deliberate failure cases.

### M9 — Explainable simulation reports

- [ ] Produce Markdown report and machine-readable findings from evidence.
- [ ] Plot physical histories alongside conservation and solver diagnostics.
- [ ] Explain timestep, theta, grid, and boundary-condition effects.
- [ ] Trace every numerical claim to a definition, dataset, and method.
- [ ] State error-estimate assumptions and limitations; do not claim bounds
      unless mathematically established.
- [ ] Provide examples of trustworthy, qualified, and failed calculations.

### M10 — Interactive numerical laboratory

- [ ] Compare synchronized physical and numerical time histories.
- [ ] Explore timestep, grid resolution, and theta interactively.
- [ ] Show linearized amplification predictions beside measured responses.
- [ ] Make explanations and evidence accessible without modifying the kernel.
- [ ] Support reproducible export of configurations and results.

### M11 — Modular physical models (original M7 scope)

- [ ] Separate EOS interface.
- [ ] Separate friction interface.
- [ ] Separate momentum-model interface.
- [ ] Separate boundary-condition interface.
- [x] Add and verify Swamee-Jain friction model.
- [ ] Add Peng-Robinson real-gas EOS closure.
- [ ] Add full momentum equation with convective acceleration.
- [ ] Enforce turbulent friction applicability and specify laminar/transition.
- [ ] Preserve constant-Z/constant-friction analytic validation oracles.
- [ ] Extend assurance evidence to every new closure.

### M12 — C++ parity and performance

- [ ] Reproduce the validated F77 algorithm in straightforward C++.
- [ ] Match state ordering and arithmetic where practical.
- [ ] Establish physical, numerical, and optimization-level parity.
- [ ] Benchmark F77 versus C++; investigate layout, aliasing, vectorization.
- [ ] Ensure both implementations emit compatible assurance evidence.

### M13 — Accelerator and surrogate research

- [ ] Identify GPU-suitable kernels without changing reference mathematics.
- [ ] Implement and validate an accelerator path against the reference.
- [ ] Benchmark CPU/GPU crossover behavior.
- [ ] Define transient quantities of interest q(xi) and sensitivities.
- [ ] Apply independently validated TDAR to pipe parameter space.
- [ ] Investigate higher-dimensional TDAR behavior.
- [ ] Preserve numerical assurance and provenance across backends.

## Future research

- [ ] Alternative spatial discretizations, including linear Galerkin.
- [ ] Alternative time integrators and adaptive timestep control.
- [ ] ThermoGPU EOS integration.
- [ ] Implicit differentiation through converged residual equations.
- [ ] Adjoint sensitivities for many-parameter problems.
- [ ] Variational/adjoint state estimation.
- [ ] Extension from a single pipe to network problems.

## Explainability and numerical assurance policy

Every simulation should return its result **and** sufficient evidence to
assess that result. Evidence generation is independent of presentation.
The numerical kernel is not responsible for prose or visualizations.

- **Proven:** follows from stated mathematics (e.g. discrete mass identity).
- **Verified:** demonstrated by tests for specified cases and tolerances.
- **Estimated:** numerical uncertainty inferred under stated assumptions.
- **Physically validated:** compared with independent physical measurements.
- **Unassessed:** evidence missing or outside a diagnostic's valid domain.

No assurance label shall imply physical validation. Every diagnostic shall
identify its mathematical definition, units, scope, threshold, evidence,
and failure behavior. Missing or invalid evidence must remain visible.

Proposed versioned run layout (subject to schema design in M8):

```text
simulation/
  results/       pressure.csv, mass_flow.csv, linepack.csv
  diagnostics/   conservation.csv, nonlinear_solver.csv,
                 discretization.csv, physical_validity.csv
  assurance/     findings.json, summary.md
  provenance/    configuration.json, numerical_methods.json
```

The M6 convergence study verifies temporal/spatial accuracy for selected
smooth cases. M7a measures errors under step, pulse, and smooth forcing;
it does **not** by itself establish damping, phase accuracy, or wave speed.
M7b explicitly addresses those remaining questions.

## Model configuration

### Default engineering models

The default physical configuration is:

- **Compressibility:** Peng-Robinson equation of state.
- **Friction:** Swamee-Jain Darcy friction factor.

The numerical solver shall access thermodynamic properties through an
EOS interface providing, at minimum,

\[
\rho(p,T,\mathbf{z})
\]

and

\[
\left(\frac{\partial \rho}{\partial p}\right)_{T,\mathbf{z}}.
\]

The pipe discretization shall not depend directly on the internal form
of the EOS.

For the Swamee-Jain model,

\[
Re = \frac{|\dot m|D}{A\mu},
\]

and

\[
f_D =
\frac{0.25}
{\left[
\log_{10}\left(
\frac{\epsilon}{3.7D}
+\frac{5.74}{Re^{0.9}}
\right)
\right]^2}.
\]

Swamee-Jain is the default engineering friction model. Its turbulent
domain of applicability shall be enforced explicitly. Laminar and
transition-flow behavior will be specified separately rather than
silently extrapolating the turbulent correlation.

### Validation models

The following simplified models are retained permanently as validation
oracles:

- constant compressibility factor \(Z\);
- constant Darcy friction factor \(f_D\).

For constant \(Z\),

\[
\rho = \frac{p}{Z R_s T},
\qquad
\frac{\partial \rho}{\partial p}
= \frac{1}{Z R_s T}.
\]

The combination of constant \(Z\) and constant \(f_D\) preserves the
analytic steady-state solution

\[
p^2(x)
=
p_{\rm in}^2
-
\frac{f_D Z R_s T}{D A^2}\dot m^2 x,
\]

and therefore remains the primary analytic regression case even after
the default engineering models become more sophisticated.

### Planned EOS implementations

- [x] Constant-Z validation EOS.
- [ ] Ideal-gas debugging EOS.
- [ ] Peng-Robinson engineering EOS.
- [ ] EOS pressure derivative verification.
- [ ] ThermoGPU/Peng-Robinson integration.

### Planned friction implementations

- [x] Constant-\(f_D\) validation model.
- [x] Swamee-Jain engineering model.
- [x] Independent Swamee-Jain reference-value tests.
- [x] Analytic \(df_D/d\dot m\).
- [x] Finite-difference verification of \(df_D/d\dot m\).
- [x] Analytic friction-source Jacobian.
- [x] Finite-difference verification of friction-source Jacobian.
- [ ] Explicit laminar-flow model.
- [ ] Explicit transition-flow policy.
