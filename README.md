# Isothermal Transient Pipe Simulator

Independent research implementation of a modular, one-dimensional isothermal transient gas-pipe simulator.

The first reference implementation will use an F77-compatible computational style. A later C++ implementation will reproduce the same numerical algorithm to study numerical and performance parity.

## V1 numerical model

- single straight, constant-diameter pipe
- fixed gas temperature and composition
- reduced isothermal momentum equation
- modular EOS and friction closures
- staggered, centered, conservative spatial discretization
- implicit theta-method time integration
- damped Newton nonlinear solve
- inlet-pressure / outlet-mass-flow reference boundary condition

The architecture is intentionally modular so that the reduced momentum model can later be replaced by the full momentum equation, and so EOS, friction, spatial discretization, time integration, boundary conditions, and solvers can be changed independently where mathematically possible.

## Repository layout

- `docs/model.tex` - governing equations, discretization, residuals, validation plan
- `docs/provenance.md` - independent-development and public-source record
- `docs/references.bib` - bibliography
- `src/fortran77/` - reference CPU implementation (planned)
- `src/cpp/` - later C++ parity implementation
- `tests/` - regression, conservation, convergence, and parity tests
- `cases/` - reproducible validation/transient cases
- `benchmarks/` - performance results and provenance

## Build the document

```bash
make docs
```

The generated PDF is `docs/model.pdf`.

## Development rule

The implementation is to be derived independently from public governing equations and public numerical literature. Prior proprietary experience may identify questions worth investigating, but is not a source for equations, algorithms, implementation details, heuristics, or test cases.
