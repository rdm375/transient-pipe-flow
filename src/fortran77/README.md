# Fortran 77-style reference implementation

Planned reference CPU solver. The computational kernel will use an F77-compatible procedural/data-layout style where practical, compiled with a modern Fortran compiler.

Planned modules/subroutines: initialization, EOS, friction, momentum contribution, spatial operators, residual assembly, Jacobian assembly, linear solve, damped Newton step, timestep advance, boundary histories, diagnostics, and output.
