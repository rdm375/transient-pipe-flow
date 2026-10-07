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
- [x] Add and verify Swamee-Jain friction model.
- [ ] Add Peng-Robinson real-gas EOS closure.
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
