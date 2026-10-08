# Transient Pipe Solver — Interface Contract

Status: M9E.3–4 checked-interface contract
Scope: Existing Fortran solver interfaces
Numerical formulation: Unchanged

## 1. Public solver entry points

The proposed public Fortran interface consists of:

- `TRANSIENT_STEP`
- `INTEGRATE_TRANSIENT`
- `INTEGRATE_TRANSIENT_ADAPTIVE`
- `INTEGRATE_TRANSIENT_ADAPTIVE_SCHEDULE`

All other routines are currently treated as implementation
details, although Fortran does not enforce this distinction.

## 2. Spatial state

For N pipe intervals:

- Pressure: PO(0:N)
- Internal mass flow: MO(0:N-1)
- Inlet mass flow: MINO

Valid N:

    2 <= N <= 100

Pressure is measured in Pa.

Mass flow is measured in kg/s.

## 3. Packed timestep solution

TRANSIENT_STEP returns a packed state U(1:2*N+2).

    U(1)       = inlet pressure
    U(2)       = inlet mass flow
    U(2*i+3)   = internal mass flow at face i
    U(2*i+4)   = pressure at node i+1

for i = 0,...,N-1.

The final pressure is U(2*N+2).

## 4. Time integration

The solver uses the theta method.

THETA is the time-weighting parameter.

The implementation accepts:

    0 <= THETA <= 1

The numerical verification suite includes:

    THETA = 0.5
    THETA = 0.65
    THETA = 1.0

These values are verified test configurations, not a restriction
on the public interface.

## 5. Required array capacities

TRANSIENT_STEP:

    PO: N+1
    MO: N
    U:  2*N+2

INTEGRATE_TRANSIENT:

    PO:    N+1
    MO:    N
    TIME:  NSTEPS+1
    HMIN:  NSTEPS+1
    HOUT:  NSTEPS+1
    HPOUT: NSTEPS+1
    HLINE: NSTEPS+1
    HNIT:  NSTEPS+1

Adaptive integrators:

    PO:    N+1
    MO:    N

    TIME, DTREC, HMIN, HOUT, HPOUT, HLINE,
    HNIT, HETA, HBAL, HCFL:

        MAXREC+1 elements each

Adaptive schedule integrator:

    TSCHED, PSCHED, QSCHED:

        NSCHED elements each

    ATOL, RTOLY:

        4 elements each

The existing assumed-size Fortran interfaces cannot verify
the actual allocation sizes supplied by callers.

## 6. Status codes

    0  Success
    1  Newton iteration limit
    2  Linear solver failure
    3  Invalid input
    4  Newton line-search failure
    5  Minimum timestep reached
    6  Adaptive rejection limit exceeded
    7  Output history capacity exhausted
    8  Invalid or excessively large adaptive indicator

Not every entry point returns every status code.

Adaptive integration may convert repeated unsuccessful
timestep attempts into status 6.

## 7. Transactional state publication

NEWTON_SOLVE:

    Copies its working state to U only after convergence.

TRANSIENT_STEP:

    Copies the converged solution to U only on success.

INTEGRATE_TRANSIENT:

    Commits physical state only after successful timesteps.

Adaptive integrators:

    Commit physical state and histories only after accepted
    timesteps.

On integration failure, the physical state represents the
last successfully accepted timestep.

## 8. History indexing

For integration histories:

    index 0 = initial condition
    index k = accepted timestep k

Adaptive integration:

    NACC is the number of accepted timesteps.

    Valid history indices after initialized execution:

        0:NACC

Fixed-step integration:

    No accepted-step count is currently returned.

    After partial failure, the caller cannot determine the
    valid history extent from INFO alone.

## 9. Validation behavior

PIPE_VALID checks:

- N in the supported range
- finite initial pressure and mass flow
- strictly positive initial pressures
- finite physical parameters
- positive DX, DT, D, A, T, Z, RS, MU
- nonnegative EPS
- positive RTOL
- nonnegative STOL
- nonnegative MAXIT
- THETA in [0,1]

Additional validation is performed by individual entry points.

Validation failure returns INFO=3.

Some diagnostic output arguments are not initialized before
early validation failures. Callers must not read undefined
outputs.

## 10. Adaptive timestep control

The adaptive integrators use a history/curvature indicator.

HETA is not a certified local truncation error estimate.

CFLMAX is an optional acoustic accuracy cap, not an implicit
stability requirement.

BFRAC optionally limits timestep size relative to boundary
forcing intervals.

The adaptive routines preserve accepted state on rejection.

## 11. Boundary forcing

Fixed-step and ramp-driven adaptive integration:

    Constant prescribed inlet pressure.
    Outlet mass flow follows a linear ramp and hold.

Adaptive schedule integration:

    Piecewise-linear prescribed inlet pressure.
    Piecewise-linear prescribed outlet mass flow.

Schedule times must be strictly increasing.

The first schedule time must be approximately zero.

The schedule must cover the requested duration.

## 12. Known interface limitations

1. Assumed-size arrays prevent independent capacity checks.

2. Fixed-step integration does not return an accepted-step
   count after partial failure.

3. Early validation failures may leave diagnostic counters
   and histories undefined.

4. Public and internal routines are not separated by a
   language-enforced interface boundary.

5. C++ interoperability is not yet defined.

## 13. Deferred numerical-behavior question

The ramp-driven adaptive integrator initializes:

    PREV(1)=PIN

The schedule-driven adaptive integrator initializes:

    PREV(1)=PO(0)

When PIN differs from PO(0), the ramp-driven initialization
includes the initial inlet-pressure discrepancy in its first
stored pressure slope.

This may affect subsequent adaptive timestep selection.

No numerical change is authorized by this document.

## 14. M9E constraints

Interface hardening must not silently change:

- the governing equations
- theta-method discretization
- staggered spatial discretization
- boundary half-cell equations
- friction or EOS calculations
- convergence tolerances
- accepted-step numerical results

Any necessary numerical change requires an independent
regression test and explicit review.

## 15. Checked-interface failure progress (M9E.3)

The `pipe_solver_api` module provides `pipe_step_checked`,
`pipe_integrate_checked`, `pipe_adaptive_checked`, and
`pipe_adaptive_schedule_checked`. Its assumed-shape arrays are
capacity-checked before dispatch to the legacy numerical kernels.

On wrapper-level capacity rejection (`INFO=3`), the physical state
and history arrays are unchanged. Adaptive checked wrappers also
return zero `NACC`, `NREJ`, `NTOTAL`, and `NNEWTON`.

After initialized adaptive execution, the physical state is the
last accepted state, even if `INFO` indicates a subsequent failure.
For the legacy arrays with lower bound zero, valid history records
are `0:NACC`. For ordinary 1-based arrays passed to the checked
module, these are elements `1:NACC+1`. Other records must not be
interpreted as results. `NTOTAL=NACC+NREJ` counts attempted steps;
`NNEWTON` includes Newton work from unsuccessful attempts.

Status 7 indicates that history capacity was exhausted after
accepted records were committed. Status 6 indicates that the
adaptive rejection limit was exceeded; an unsuccessful attempt
does not commit the physical state. Early numerical-input
validation inside the legacy kernel may leave histories and
diagnostic outputs undefined; the wrapper's zero counters do not
constitute a guarantee about the kernel's early-failure outputs.

Fixed-step integration has no accepted-step count in its existing
interface. Do not infer a valid partial-history extent from `INFO`
alone. Adding such a count requires a separate API revision.

## 16. Build isolation (M9E.4)

Debug regression binaries use shared checked-build object files.
Optimized characterization programs remain separate compilation
paths and must not link debug objects. The M9C.3 failure-injection
regression substitutes a test `TRANSIENT_STEP` and must never link
the production timestep object.
