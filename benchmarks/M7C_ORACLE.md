# M7c experiment 2: exact-defect oracle and propagation separation

This is a **diagnostic of the linearized semidiscrete system**, not a
production global-error estimator. It uses the M7b Jacobian export and
`scipy.linalg.expm` to distinguish local-defect approximation from error
transport. It performs no step doubling and changes no Fortran code.

For `y_num_next = G y_num` and `y_exact_next = E y_exact`, with
`G = (I - theta*h*J)^-1 (I + (1-theta)*h*J)` and `E = exp(h*J)`, define
`error = y_num - y_exact` and `defect = (G-E)*y_num`.
The **exact** identity is

    error_next = E*error + defect.

Notice that the propagation matrix is **E**, not G, when the defect is
evaluated at the numerical state. Using G instead is an approximation;
using G with the exact defect is not an exact oracle identity.

The script compares five cases:

* `oracle-exact`: exact defect and E propagation (identity check).
* `curvature-exact`: `(theta-.5)*h^2*J^2*y_num`, E propagation.
* `resolvent-exact`: implicit-resolvent filtered curvature, E propagation.
* `curvature-theta`: curvature defect, G propagation.
* `resolvent-theta`: filtered curvature defect, G propagation; matches the
  structure of the first M7c experiment.

Each row includes the local defect's relative error, final vector relative
error, signed alignment, and effectivity. The oracle identity is checked
against a `1e-7` relative threshold (actual results should be much better).

Run from the repository root:

    OPENBLAS_NUM_THREADS=1 python benchmarks/global_error_oracle.py

Outputs: `benchmarks/results-m7c/oracle_history.csv` and
`benchmarks/results-m7c/oracle_summary.csv`.

A zero final error is not evidence of an accurate estimator: ratios can be
ill-conditioned when the actual final error is small. Inspect absolute
errors and histories as well as effectivity.
