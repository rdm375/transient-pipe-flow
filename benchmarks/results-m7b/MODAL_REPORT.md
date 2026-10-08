# M7b — Linearized modal verification

The operator is constructed directly from the production F77 analytical
Jacobian evaluated at the analytical steady state, with fixed inlet
pressure and prescribed outlet flow. The inlet-flow algebraic constraint
is eliminated. The remaining 2N variables are face flows and interior
nodal pressures. Thus the modal results use the **actual discrete**
staggered geometry, including the outlet half-cell.

For each continuous eigenvalue λ of the semidiscrete operator,
the theta integrator predicts g=(1+(1−θ)hλ)/(1−θhλ).
The decay rate is −log|g|/h and the sampled angular frequency
is arg(g)/h. Frequencies are aliased modulo 2π/h.

| N | dt (s) | θ | max |g| | min |g| |
|---:|---:|---:|---:|---:|
| 20 | 15 | 0.50 | 0.976852 | 0.735347 |
| 20 | 15 | 0.65 | 0.976932 | 0.676016 |
| 20 | 15 | 1.00 | 0.977117 | 0.408539 |
| 20 | 60 | 0.50 | 0.966866 | 0.242214 |
| 20 | 60 | 0.65 | 0.911701 | 0.319558 |
| 20 | 60 | 1.00 | 0.914348 | 0.114039 |
| 20 | 120 | 0.50 | 0.982641 | 0.099123 |
| 20 | 120 | 0.65 | 0.832987 | 0.056441 |
| 20 | 120 | 1.00 | 0.842211 | 0.057550 |
| 40 | 15 | 0.50 | 0.976842 | 0.735355 |
| 40 | 15 | 0.65 | 0.976923 | 0.590246 |
| 40 | 15 | 1.00 | 0.977108 | 0.223634 |
| 40 | 60 | 0.50 | 0.991295 | 0.242231 |
| 40 | 60 | 0.65 | 0.911666 | 0.319572 |
| 40 | 60 | 1.00 | 0.914315 | 0.057643 |
| 40 | 120 | 0.50 | 0.995594 | 0.099105 |
| 40 | 120 | 0.65 | 0.832924 | 0.056454 |
| 40 | 120 | 1.00 | 0.842154 | 0.028890 |

## Interpretation and boundaries

- Eigenmode damping is a **linearized prediction**, not a fitted nonlinear response.
- The sensor histories include nonlinear friction, finite disturbances, and boundary forcing.
- A slope reversal does not identify a particular eigenmode or prove instability.
- Step/pulse endpoint quadrature can introduce forcing errors distinct from modal damping.
- This calculation does not establish continuum wave speed, a nonlinear error bound, or physical validation.
- A direct small-perturbation experiment is still required to establish modal agreement with nonlinear histories.
- `modal_metrics.csv` contains every complex eigenmode at each N, dt and theta.
