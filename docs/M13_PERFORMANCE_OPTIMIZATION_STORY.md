# M13 — CPU Performance Optimization Story

**Project:** Isothermal transient single-pipe simulator  
**Status:** COMPLETE  
**Final milestone:** `m13-final`  
**Platform:** Dell Precision 7710, Debian Linux  
**Compiler:** GCC/GFortran 15.2  
**Optimization:** `-O3 -march=native -flto`

## 1. Objective

M13 investigated whether an independently developed C++20
implementation of the transient pipeline simulator could match
the performance of the original Fortran77 implementation.

The acceptance criterion was:

    median(C++ CPU time) / median(Fortran CPU time) <= 1.05

Numerical agreement was a prerequisite. Performance improvements
were not permitted to conceal unexplained numerical discrepancies
by widening tolerances.

The reference workload was:

| Parameter | Value |
|---|---:|
| Pipe length | 100 km |
| Spatial cells | 100 |
| Spatial resolution | 1 km |
| Timestep | 15 s |
| Timesteps | 8,000 |
| Simulated duration | 120,000 s |
| Total Newton iterations | 1,049 |

## 2. Initial performance investigation

Early measurements suggested that C++ was faster than Fortran.

However, these measurements used executables with different
optimization histories.

Controlled rebuilding revealed that compiler configuration
materially affected the comparison.

Without link-time optimization, C++ was faster than the
rebuilt Fortran implementation.

With `-O3 -march=native -flto`, Fortran became faster.

The original matched-LTO comparison measured:

| Implementation | Median CPU time |
|---|---:|
| Fortran77 | 131.778 ms |
| C++20 | 142.291 ms |

The C++/Fortran ratio was 1.07978.

This exceeded the acceptance threshold.

The result demonstrated that conclusions about language
performance could not be separated from build configuration.

## 3. Investigation of the performance difference

Callgrind profiling identified approximately:

- 660.9 million instruction references for Fortran.
- 742.4 million instruction references for C++.

C++ executed approximately 12.34% more instructions.

The major mathematical-library routines had identical
exclusive instruction counts in the measured profiles.

Consequently, the difference could not be explained by
additional execution inside those particular libm routines.

The investigation examined:

- Banded matrix indexing.
- Jacobian assembly.
- Newton iteration overhead.
- Residual evaluation.
- Workspace allocation and reuse.
- Repeated physical-model calculations.
- Compiler inlining and specialization.

These investigations motivated several targeted optimizations.

## 4. C++ M13.2–M13.4

Three important C++ improvements were introduced.

### M13.2b — Direct-indexed banded solver

The banded linear solver was specialized for the
production Jacobian bandwidth.

This reduced indexing overhead while preserving
the numerical elimination algorithm.

### M13.3 — Old-time friction-source caching

The old-time friction contribution is invariant
during Newton iterations within a timestep.

Caching this contribution eliminated repeated
evaluations during residual assembly.

### M13.4 — Workspace reuse

The initial Newton guess was moved into
preallocated workspace storage.

This reduced repeated memory-management overhead.

After these optimizations, a controlled comparison measured:

| Implementation | Median CPU time |
|---|---:|
| Fortran LTO | 131.974 ms |
| C++ M13.4 | 137.375 ms |

The C++/Fortran ratio was 1.040921.

The original acceptance criterion was met.

However, this comparison used the existing Fortran implementation.
It did not establish parity against an independently optimized
Fortran baseline.

## 5. Symmetric Fortran optimization

The next objective was to apply comparable optimization
effort to Fortran.

### M13.5 — Old-time friction-source caching

The Fortran residual implementation was modified to cache
old-time friction contributions.

The resulting end-to-end performance improvement was small
and comparable to measurement variability.

### M13.6 — Specialized banded solver

A direct-indexed Fortran banded solver was developed
and tested separately.

The specialized solver demonstrated approximately
1.304x speedup in its isolated microbenchmark.

However, the full simulation improved by only
approximately 0.7%, a noise-level result.

The experimental solver was not integrated into
the production implementation.

This demonstrated that microbenchmark improvements
do not necessarily translate into application-level gains.

## 6. M13.7 — The major Fortran breakthrough

Profiling identified repeated Swamee–Jain friction-factor
evaluations as a significant cost.

The key observation concerned the Newton initial guess.

At the beginning of a timestep, the initial mass-flow
guess equals the previous timestep's mass flow.

Therefore, the Reynolds number and friction factor
are unchanged at the initial guess.

The old-time friction factor can be reused.

However, the old-time friction source cannot generally
be reused because the new-time face density may differ.

The optimization therefore:

1. Calculates the old-time friction factor once.
2. Stores the factor in the Newton workspace.
3. Reuses it when the initial mass flow matches
   the old-time mass flow.
4. Recomputes the friction source using the current density.
5. Retains normal friction evaluation when flows differ.

The resulting controlled benchmark measured:

| Fortran implementation | Median CPU time |
|---|---:|
| Original LTO | 130.759 ms |
| M13.7 optimized | 90.662 ms |

The optimization reduced CPU time by approximately 30.7%.

The full 8,000-step simulation retained its validated
numerical results.

## 7. M13.8 — Equivalent C++ optimization

The same optimization was implemented in C++20.

The Newton workspace gained storage for old-time
friction factors.

The initial residual evaluation uses a specialized path
that reuses these factors when the mass flows match.

The existing line-search residual evaluation remains
unchanged.

Numerical validation included:

- Successful C++ LTO compilation.
- Five passing CTest regression tests.
- 606 exact residual comparisons.
- Exact agreement of seven shared reported numerical
  quantities in the full 8,000-step simulation.

The final 60-round alternating CPU benchmark measured:

| Implementation | Median CPU time |
|---|---:|
| Fortran M13.7 | 90.998 ms |
| C++ M13.4 | 137.518 ms |
| C++ M13.8 | 95.394 ms |

C++ M13.8 achieved:

- 30.632% lower CPU time than C++ M13.4.
- 1.4416x speedup over C++ M13.4.
- C++/Fortran CPU-time ratio of 1.048309.

The final performance criterion was satisfied.

The margin was narrow, and the result is specific
to the measured workload and hardware.

## 8. Engineering lessons

### 8.1 Compare equivalent build configurations

Compiler optimization settings can reverse
the apparent performance ranking.

### 8.2 Optimize both implementations

Matching an unoptimized reference is not sufficient
to establish symmetric performance parity.

### 8.3 Measure application-level effects

The M13.6 banded solver microbenchmark improved
substantially, but its effect on total runtime was small.

### 8.4 Exploit mathematical invariants

The largest performance improvement came from
recognizing that the initial Newton mass flow
equals the old-time mass flow.

This eliminated redundant physical-model work
without changing the governing equations.

### 8.5 Preserve numerical correctness

The optimization reused friction factors rather
than incorrectly reusing complete friction sources.

The distinction matters because density can change
even when mass flow remains unchanged.

### 8.6 Preserve reproducibility

The final benchmark recorded executable checksums,
workload parameters, numerical outputs, and
external process CPU-time measurements.

## 9. Final result

| Metric | Fortran77 | C++20 |
|---|---:|---:|
| Median CPU time | 90.998 ms | 95.394 ms |
| CPU time per timestep | 11.375 us | 11.924 us |
| Timesteps | 8,000 | 8,000 |
| Newton iterations | 1,049 | 1,049 |

Final performance ratio:

    C++ / Fortran = 1.048309

Acceptance threshold:

    C++ / Fortran <= 1.05

**M13 PERFORMANCE ACCEPTANCE: PASS**

## 10. Future work

M13 establishes validated CPU reference implementations.

Further optimization is deliberately deferred.

The remaining CPU performance ceiling will be investigated
separately under M15.

GPU execution models and accelerator feasibility
will be investigated under M14.

The M13 implementations and numerical results remain
the reference baseline for both investigations.

**M13 CLOSED.**
