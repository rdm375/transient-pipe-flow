# M7c exploratory experiment: Jacobian-based global error transport

This experiment is **linearized** and uses the M7b matrix exported from
production F77 Jacobian code. It is not a nonlinear global-error estimate.
It does not use step doubling.

Run with:

```sh
OPENBLAS_NUM_THREADS=1 python benchmarks/global_error_transport.py
```

Requires NumPy and SciPy. The input is
`benchmarks/results-m7b/modal_matrix.csv`. Output is
`benchmarks/results-m7c/transport_history.csv` and
`transport_summary.csv`.

For `u' = J u`, the theta map is

`G = (I - theta*h*J)^(-1) (I + (1-theta)*h*J)`.

The experiment predicts numerical-minus-exact global error using

`e_hat[n+1] = G*e_hat[n] + tau_hat[n]`,

`tau_hat[n] = (theta-1/2)*h^2*(I-theta*h*J)^(-1)*J^2*u[n]`.

The resolvent modification is a **heuristic** to improve behavior for stiff
modes, not an exact local-defect identity or a proven error bound.
The reference is `expm(h*J)` applied repeatedly to the same initial state.
The initial condition is a unit perturbation in a central pressure coordinate;
the operator is in the M7b flow/pressure interleaved ordering. Euclidean norms
are **unweighted**, so they are not engineering tolerance norms.

The estimator's effectivity varies substantially, especially at coarse
steps for theta=0.65. This negative result motivates reconstruction-based
residual transport rather than further tuning an unverified curvature factor.

Next: construct a continuous reconstruction and evaluate its semidiscrete
defect at quadrature nodes. Transport that defect through the tangent system.
Then compare signed vector errors and physical pressure/flow output errors
against independent nonlinear reference trajectories.
