# Future Numerical Accuracy and Robustness Refinements

**Status:** Deferred
**Origin:** M9 Newton convergence investigation, October 2026
**Priority:** Revisit when justified by accuracy, robustness, or performance requirements.

## 1. Background

During M9 development, a fixed-timestep simulation encountered Newton
line-search stagnation.

Test configuration:

- Pipe length: 100 km
- Pipe diameter: 1 m
- Spatial cells: N = 20
- Timestep: 1.875 s
- Simulation duration: 3600 s
- Theta-method parameter: 0.65
- Constant-Z equation of state

The original solver failed at timestep 1592, near t = 2985 s.

Diagnostics indicated:

- Maximum continuity residual: approximately 1.34e-11 kg/s
- Newton residual tolerance: 1e-11
- Pressure corrections approaching floating-point resolution
- Strict line search unable to obtain further residual reduction

The evidence is consistent with floating-point stagnation rather than
a physically significant failure of the transient model.

## 2. Experimental Guarded Newton Convergence

An experimental convergence rule was implemented without modifying
the production solver.

The rule requires:

1. The original residual to remain within 10 times RTOL.
2. Independently scaled continuity and momentum residuals to satisfy
   prescribed accuracy limits.
3. Pressure and mass-flow corrections to satisfy independently scaled
   convergence criteria.
4. The pressure boundary residual to satisfy its own criterion.

Results:

| N | dt (s) | Scaled tolerance | Result | Guarded acceptances |
|---|--------|------------------|--------|---------------------|
| 20 | 1.875 | 1e-13 | Failed | 0 |
| 20 | 3.750 | 1e-13 | Passed | 0 |
| 20 | 1.875 | 1e-12 | Passed | 8 |
| 20 | 3.750 | 1e-12 | Passed | 0 |
| 20 | 1.875 | 1e-11 | Passed | 8 |
| 20 | 3.750 | 1e-11 | Passed | 0 |
| 40 | 1.875 | 1e-12 | Passed | 0 |
| 40 | 3.750 | 1e-12 | Passed | 0 |
| 80 | 1.875 | 1e-12 | Passed | 0 |
| 80 | 3.750 | 1e-12 | Passed | 0 |

For N = 20 and dt = 1.875 s, the guarded experiment completed
all 1920 timesteps.

The final cumulative numerical mass-conservation defect was
approximately -5.35e-10 kg.

These results establish feasibility for the tested scenario, not
general reliability across operating conditions.

The guarded convergence rule has not been adopted as the production
default.

## 3. Potential Floating-Point Refinements

### 3.1 Density-storage increments

For the constant-Z EOS,

    rho(p) = p / (Z Rs T)

evaluate the density change as

    delta_rho = (p_new - p_old) / (Z Rs T)

rather than subtracting independently evaluated densities.

This is algebraically equivalent in exact arithmetic and may reduce
rounding error in the continuity storage term.

Apply consistently to interior cells and boundary half-cells.

### 3.2 Theta-weighted mass-flow differences

Evaluate

    theta * (m_right_new - m_left_new)
      + (1 - theta) * (m_right_old - m_left_old)

rather than independently computing theta-weighted flows and then
subtracting them.

The benefit must be measured against a higher-precision reference.

### 3.3 Direct inventory increments

For the constant-Z EOS, compute the inventory change directly from
pressure increments:

    delta_M = A dx / (Z Rs T)
              * sum_i w_i * (p_i_new - p_i_old)

where w_i = 1 for interior nodes and 1/2 at the endpoints.

Compare this against

    M_new - M_old

to distinguish rounding in inventory diagnostics from numerical
conservation defects.

Preserve the distinction between:

- Physical mass imbalance: inlet flow minus outlet flow.
- Inventory change rate.
- Numerical mass-conservation defect.

These quantities must not be conflated.

### 3.4 Momentum residual evaluation

Investigate rounding and cancellation in combinations of:

- Pressure gradients.
- Transient mass-flow terms.
- Friction contributions.

Evaluate alternative expression groupings against a higher-precision
reference before changing the production implementation.

### 3.5 Jacobian scaling and linear solver

Investigate:

- Row and column equilibration.
- Scaled pivot tests.
- Linear-system backward errors.
- Conditioning across spatial resolutions and timesteps.

Do not modify pivot tolerances without numerical evidence.

### 3.6 Newton state updates

Investigate:

- Componentwise convergence criteria.
- Floating-point update-resolution detection.
- Guarded stagnation acceptance.
- Pressure and mass-flow scaling.
- Increment-based state variables.

Any convergence safeguard must reject genuinely unconverged states.

### 3.7 Adaptive timestep error estimates

Investigate roundoff sensitivity in:

- Differences between successive temporal slopes.
- Inventory changes.
- Small physical mass imbalances.
- Near-steady-state error estimates.

A roundoff-dominated estimator may cause unnecessary timestep
rejection or inappropriate timestep selection.

### 3.8 Compensated and extended-precision arithmetic

Potential techniques include:

- Compensated summation for inventory accumulation.
- Higher-precision reference residual evaluation.
- Selective extended precision for sensitive calculations.
- Mixed-precision strategies for future accelerator implementations.

These are research options, not current requirements.

## 4. Proposed Evaluation Method

Before adopting any reformulation:

1. Preserve the original implementation as a control.
2. Evaluate old and new expressions at identical stored states.
3. Compare both against higher-precision reference evaluations.
4. Verify algebraic equivalence of the discrete equations.
5. Verify consistency with the analytic Jacobian.
6. Run fixed-timestep and adaptive regression tests.
7. Compare pressure, mass flow, inventory, and physical imbalance.
8. Measure local and integrated numerical conservation defects.
9. Evaluate Newton convergence and runtime.
10. Reject unexplained numerical discrepancies rather than widening
    tolerances.

A smaller double-precision residual alone does not prove that the
physical solution is more accurate.

## 5. Suggested Priority

### Tier 1: Low-complexity numerical improvements

- Density-storage increments.
- Theta-weighted flow differences.
- Direct inventory increments.

### Tier 2: Solver robustness

- Momentum residual evaluation.
- Jacobian equilibration.
- Linear-solver backward-error diagnostics.
- Adaptive estimator roundoff analysis.

### Tier 3: Research investigations

- Increment-based Newton formulation.
- Compensated arithmetic.
- Extended and mixed precision.

## 6. Decision

These investigations are deferred.

The current simulator already demonstrates strong numerical
mass-conservation performance in the tested cases.

Future development should prioritize the main simulator roadmap.

Revisit these refinements if:

- Newton stagnation prevents required simulations.
- Smaller timesteps expose a residual accuracy floor.
- Conservation accuracy becomes inadequate.
- Adaptive timestep control becomes roundoff-limited.
- New EOS models or accelerator backends expose additional
  floating-point sensitivity.

No production algorithm or tolerance change is authorized by
this document.
