# Isothermal Pipe Simulator Milestones

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

- [ ] Refactor residual evaluation into reusable F77-style routines.
- [ ] Define the complete transient unknown-vector ordering.
- [ ] Implement complete theta-method transient residual.
- [ ] Derive analytic Jacobian entries.
- [ ] Implement analytic Jacobian assembly.
- [ ] Implement independent finite-difference Jacobian.
- [ ] Compare analytic and finite-difference Jacobians.
- [ ] Test Jacobian at non-equilibrium states.
- [ ] Verify expected Jacobian sparsity/locality.

### M3 — Newton nonlinear solve

- [ ] Implement linear solve for the reference Jacobian.
- [ ] Implement Newton iteration.
- [ ] Implement residual-based convergence criteria.
- [ ] Implement damped Newton/backtracking.
- [ ] Verify recovery of analytic steady state from perturbed guesses.
- [ ] Verify failure reporting for nonconvergent cases.

### M4 — Single time step

- [ ] Implement one theta-method transient step.
- [ ] Preserve an exact steady state for one step.
- [ ] Verify per-step discrete global mass balance.
- [ ] Verify boundary histories at t^n and t^(n+1).
- [ ] Verify inlet boundary flow diagnostic.

### M5 — Transient integration

- [ ] Implement multi-step integration.
- [ ] Implement downstream mass-flow ramp case.
- [ ] Record pressure, mass-flow, linepack, and boundary histories.
- [ ] Verify global mass conservation over complete simulation.
- [ ] Verify relaxation to the analytic steady solution.

### M6 — Numerical characterization

- [ ] Temporal refinement study.
- [ ] Spatial refinement study.
- [ ] Theta study: 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 1.00.
- [ ] Characterize numerical damping and transient accuracy.
- [ ] Establish reference cases and tolerances.

### M7 — Modular physical models

- [ ] Separate EOS interface.
- [ ] Separate friction interface.
- [ ] Separate momentum-model interface.
- [ ] Separate boundary-condition interface.
- [ ] Add Reynolds-dependent friction.
- [ ] Add real-gas EOS closure.
- [ ] Add full momentum equation with convective acceleration.

## V2 — C++ parity implementation

- [ ] Reproduce F77 numerical algorithm in straightforward C++.
- [ ] Match state ordering and arithmetic where practical.
- [ ] Establish physical parity.
- [ ] Establish numerical parity.
- [ ] Establish compiler/optimization parity.
- [ ] Benchmark F77 versus C++.
- [ ] Investigate vectorization, aliasing, layout, and compiler effects.

## V3 — Accelerator and surrogate work

- [ ] Identify GPU-suitable kernels.
- [ ] Implement accelerator path without changing reference mathematics.
- [ ] Benchmark CPU/GPU crossover behavior.
- [ ] Define transient quantities of interest q(xi).
- [ ] Compute parameter sensitivities.
- [ ] Apply TDAR to the transient-pipe parameter space.
- [ ] Investigate higher-dimensional TDAR behavior.

## Future research

- [ ] Alternative spatial discretizations, including linear Galerkin.
- [ ] Alternative time integrators.
- [ ] ThermoGPU EOS integration.
- [ ] Implicit differentiation through converged residual equations.
- [ ] Adjoint sensitivities for many-parameter problems.
- [ ] Variational/adjoint state estimation.
- [ ] Extension from a single pipe to network problems.
