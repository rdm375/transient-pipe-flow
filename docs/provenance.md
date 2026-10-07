# Provenance and Independent Development

## Policy

This project is independently developed from publicly available literature and first-principles governing equations. No proprietary source code, internal documentation, proprietary algorithms, proprietary test cases, or implementation materials from commercial pipeline simulators are to be used.

Prior industry experience may identify a physical or numerical question worth investigating. The solution used here must then be independently derived or supported by a public source and documented below.

## Initial public lineage

| Component | Public source | Use in this project |
|---|---|---|
| 1-D isothermal conservation equations | Kranjčić et al. (2024) and standard compressible-flow conservation laws | Public statement of mass/momentum equations and assumptions; equations are re-derived in `model.tex`. |
| Implicit finite differences and Newton solution | Kiuchi (1994) | Establishes public precedent for fully implicit finite-difference transient gas simulation solved with Newton-Raphson. |
| Staggered finite differences, general EOS, mass conservation | Gyrya & Zlotnik (2019) | Establishes public precedent for staggered pressure/density/flow discretization and exact discrete mass conservation. Our implicit theta formulation is derived independently. |
| Theta time integration | Kranjčić et al. (2024) | Establishes public use of a finite-difference theta-scheme for transient pipeline dynamics. Our residual equations are derived independently. |

## V1 independently selected/derived choices

- Single horizontal, constant-diameter pipe.
- Reduced momentum equation initially; full convective momentum retained as a future interchangeable model.
- Pressure/density at nodes and mass flow at staggered face locations.
- Centered differences in space.
- Arithmetic node-to-face density interpolation for V1. For constant Z and constant Darcy friction this reproduces the analytic steady pressure-squared relation at grid nodes.
- Theta is a configurable numerical parameter; no undocumented remembered default is to be adopted.
- Constant-Z and constant-friction reference case is the first analytic validation target.

## Review rule

If a design choice appears unusually similar to a remembered proprietary implementation, stop and document the concern. Search the public literature/patent record before proceeding. Either establish an independent public lineage or choose a different formulation.
