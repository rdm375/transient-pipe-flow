# M7b — Small-perturbation linearization verification

The linear model uses the production analytical Jacobian and the
outlet half-cell forcing derivative. Both models use the same
theta-method and exactly aligned smooth demand forcing.

The comparison is made at matching timesteps, so this test examines
linearization consistency, **not temporal discretization accuracy**.

| dt (s) | theta | pressure error at amplitude 1 (Pa) | pressure error at amplitude 0.125 (Pa) | observed order |
|---:|---:|---:|---:|---:|
| 15 | 0.50 | 8.51612 | 0.1333405 | 1.999 |
| 15 | 0.65 | 8.512855 | 0.1332884 | 1.999 |
| 15 | 1.00 | 8.505511 | 0.1331709 | 1.999 |
| 60 | 0.50 | 8.514111 | 0.1333092 | 1.999 |
| 60 | 0.65 | 8.501236 | 0.1331034 | 1.999 |
| 60 | 1.00 | 8.476243 | 0.1327029 | 1.999 |

## Interpretation

- Quadratic discrepancy in perturbation amplitude is expected from
  a smooth nonlinear model with a consistent first-order tangent.
- Very small perturbations eventually encounter roundoff and Newton
  convergence tolerance floors; measured orders need not be exactly 2.
- Agreement verifies the linearized forced time-step response, not
  a decomposition of finite-amplitude histories into individual modes.
- The eigenvalues of the same operator independently determine the
  theta amplification factors reported in `MODAL_REPORT.md`.
- This does not verify continuum wave speed or physical validity.
