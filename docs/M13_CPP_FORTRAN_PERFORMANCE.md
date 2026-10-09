# M13 — C++20 / Fortran77 performance parity investigation

**Project:** Isothermal transient single-pipe simulator (`isothermal-pipe-sim`)
**Status:** C++ performance target met against the existing Fortran LTO baseline; symmetric Fortran optimization and re-baselining pending
**Platform:** Dell Precision 7710, GNU GCC/GFortran 15.2, Debian Linux
**Performance acceptance criterion:** `median(C++ CPU time) / median(Fortran CPU time) <= 1.05`
**Reference workload:** 100 spatial cells, 15 s timestep, 8,000 steps (120,000 s simulated)

## Executive summary

The project ported a validated fixed-form Fortran77 transient gas-pipeline solver to native C++20, retaining the numerical formulation and testing each translated component against its Fortran counterpart. Initial measurements suggested that C++ was faster than Fortran. Controlled recompilation revealed that this result depended strongly on **link-time optimization (LTO)**: C++ led in matched non-LTO builds, whereas Fortran led when both implementations were compiled with `-O3 -march=native -flto`.

The original matched-LTO experiment found medians of **131.778 ms (Fortran)** and **142.291 ms (C++)**, a C++/Fortran ratio of **1.07978**. After targeted C++ optimizations M13.2b–M13.4, a new 40-round rotating-order external process-CPU-time experiment measured **131.974 ms (Fortran LTO)**, **142.957 ms (original C++)**, and **137.375 ms (optimized C++)**. The optimized ratio was **1.040921**, meeting the 1.05 acceptance threshold with **1.198 ms** margin. This is parity with the **existing** Fortran LTO baseline, not a claim about an independently reoptimized Fortran implementation. The ratio is workload- and machine-specific; it does not rank languages generally.

Crucially, both LTO executables produced matching reported final outlet pressure, linepack, Newton iteration count, and maximum numerical mass-conservation defect. Callgrind reported **660,855,731 instructions** for Fortran and **742,413,263** for C++, a difference of **81,557,532** (12.34% relative to Fortran). The dominant measured libm routines had **identical exclusive instruction counts**, so the excess is not explained by additional execution *inside those routines*. The historical instruction difference motivated targeted changes; the results of those changes and remaining uncertainties are documented in Sections 11–14.

## 1. Scientific and numerical context

The simulator solves isothermal transient flow in one pipe using a staggered discretization: nodal pressure/density, staggered mass flow, and a theta-method time integrator. The benchmark prescribes inlet pressure and an outlet demand ramp from 100 to 105 kg/s over 600 s. The Newton method uses a banded Jacobian and banded linear solve. The C++ port preserves the Fortran arithmetic grouping where required for numerical parity.

The model tracks linepack and a discrete conservation residual. **Physical mass-flow imbalance** (inlet minus outlet flow) is not the same quantity as **numerical mass-conservation defect** (the mismatch between inventory change and theta-weighted boundary flow imbalance). The reported `max_mass_defect` is the latter.

### Translation and validation

The port was implemented incrementally, with parity tests for physical models, residual, Jacobian, banded solver, Newton solver, transient timestep, and integration. The accumulated comparison count reached **80,533** across these test suites, including the 8,000-step integration parity case. Earlier parity runs showed zero observed relative discrepancy for the compared quantities under the tested compiler settings. This does **not** establish bitwise identity of all intermediate states for every compiler configuration.

The CMake C++ executable has no dependency on `libgfortran`. The profiling C++ LTO build passed all five CTest checks:

1. `cpp_short_transient`
2. `cpp_long_transient`
3. `cpp_zero_steps`
4. `cpp_numerical_reference_80`
5. `cpp_numerical_reference_8000`

The numerical-reference tests are tolerance-based; the observed printed equality of selected final values is a separate observation.

## 2. Benchmark methodology and build configurations

### Workload

- Spatial cells: `n=100`.
- Time step: `dt=15 s`.
- Simulation duration: `120000 s`.
- Time steps: `8000`.
- Total Newton iterations: `1049` in both implementations.
- Outlet demand ramp: 100 to 105 kg/s over 600 s.

The Fortran benchmark accepts `100 15 120000`; the C++ benchmark accepts `100 8000`. Both correspond to the same numerical case.

### Compilation

| Build | Relevant optimization flags | Notes |
|---|---|---|
| Fortran baseline | `-O3 -march=native` | Rebuilt via `make -B m10-banded-build` |
| C++ baseline | `-O3 -march=native` | CMake native build |
| Fortran LTO | `-O3 -march=native -flto` | `make -B FC=gfortran FFLAGS_OPT="-O3 -march=native -flto" build/m10_cpu_banded` |
| C++ LTO | `-O3 -march=native -flto` | Separate CMake build `build/m13-cpp-lto` |
| C++ LTO + symbols | `-O3 -march=native -flto -g` | Separate CMake build `build/m13-cpp-profile` |

**Caveat:** matching top-level optimization flags does not guarantee equivalent optimization semantics across C++ and Fortran. Language rules, source representation, aliasing, runtime checks, and GCC front ends differ. Also, the Fortran benchmark uses `CPU_TIME` internally while the C++ application reports elapsed time; the final controlled comparison instead measured both child processes externally in the same way.

### Timing procedure

Repeated runs were alternated to reduce simple ordering bias. A later experiment used Python `os.fork()`, `os.execv()`, and `os.wait4()` to obtain each child process's user-plus-system CPU time; `time.perf_counter()` simultaneously recorded external wall time. Thirty repetitions per implementation were collected after warm-up, with alternating order. The process-CPU measurement includes startup and shutdown, unlike a narrowly timed integration kernel. Results are medians; no confidence intervals or CPU-frequency/thermal controls have yet been reported.

## 3. Performance chronology

### 3.1 Initial comparison: C++ appears faster

A 20-run alternating external wall-time comparison against a preserved M11 Fortran executable yielded:

| Executable | Median | Minimum | Maximum |
|---|---:|---:|---:|
| Preserved Fortran reference | 184.519 ms | 179.998 ms | 199.594 ms |
| C++ native, non-LTO | 155.581 ms | 152.366 ms | 166.041 ms |

The ratio `184.519 / 155.581 = 1.1860` suggests C++ was approximately **1.186× as fast** for this comparison. The precise historical flags of the preserved Fortran binary could not be verified from the executable alone, so this was not the strongest controlled result.

### 3.2 Rebuilt non-LTO comparison

A subsequent 15-run, three-way comparison yielded:

| Executable | Median | Minimum | Maximum |
|---|---:|---:|---:|
| Preserved M11 Fortran | 182.750 ms | 179.958 ms | 228.787 ms |
| Rebuilt Fortran, non-LTO | 175.914 ms | 171.997 ms | 191.182 ms |
| C++ native, non-LTO | 153.892 ms | 152.324 ms | 165.414 ms |

The C++/rebuilt-Fortran time ratio was `153.892 / 175.914 = 0.8748`, or a C++ speedup of approximately **1.143×** (about **12.5% lower elapsed time**). This confirmed that the initial C++ advantage was not solely an artifact of the preserved executable.

### 3.3 Fortran LTO changes the ranking

After rebuilding Fortran with `-flto`, a 20-run alternating comparison yielded:

| Executable | Median | Minimum | Maximum |
|---|---:|---:|---:|
| Fortran LTO | 131.375 ms | 128.140 ms | 161.778 ms |
| C++ non-LTO | 153.469 ms | 151.964 ms | 213.694 ms |

The performance ranking reversed. This established that the earlier language-level conclusion was premature: **build configuration materially changed the comparison**.

### 3.4 Matched LTO builds

C++ was then rebuilt with `-flto`. A 20-run, three-way comparison yielded:

| Executable | Median | Minimum | Maximum |
|---|---:|---:|---:|
| Fortran LTO | 133.206 ms | 127.659 ms | 156.703 ms |
| C++ non-LTO | 154.998 ms | 151.190 ms | 178.824 ms |
| C++ LTO | 142.500 ms | 138.300 ms | 152.519 ms |

The matched-LTO ratio was `142.500 / 133.206 = 1.0698`: C++ was about **6.98% slower**. LTO improved both builds, but improved the Fortran executable substantially more relative to its non-LTO baseline. The exact causal transformations responsible were not yet established.

### 3.5 External CPU time confirms the remaining gap

A 30-run alternating experiment measured both external wall and process CPU time:

| Metric | Fortran LTO | C++ LTO |
|---|---:|---:|
| Median wall time | 132.316 ms | 142.673 ms |
| Minimum wall time | 129.761 ms | 138.497 ms |
| Maximum wall time | 149.300 ms | 162.529 ms |
| Median CPU time | **131.778 ms** | **142.291 ms** |
| Minimum CPU time | 129.363 ms | 138.147 ms |
| Maximum CPU time | 148.944 ms | 162.144 ms |

The CPU-time ratio was **1.07978**. Wall and CPU medians were close, indicating that scheduler waiting did not dominate the observed difference. The experiment does not isolate clock-frequency changes, memory effects, or microarchitectural causes.

**M13 target:** `1.05 * 131.778 = 138.3669 ms`; C++ must improve by approximately `142.291 - 138.3669 = 3.9241 ms` relative to this baseline.

## 4. Numerical equivalence under LTO

Both validated LTO executables completed the full benchmark. Representative outputs:

| Quantity | Fortran LTO | C++ LTO |
|---|---:|---:|
| Steps | 8000 | 8000 |
| Duration | 120000 s | 120000 s |
| Newton iterations | 1049 | 1049 |
| Final outlet pressure | 7,843,657.8677442567 Pa | 7,843,657.8677442567 Pa |
| Final linepack | 4,798,420.4378263094 kg | 4,798,420.4378263094 kg |
| Maximum mass-conservation defect | 6.2161021219253598e-08 kg/s | 6.2161021219253598e-08 kg/s |

C++ reported final inlet flow `104.99999993783898 kg/s` and final outlet flow `105 kg/s`. Their difference is approximately `-6.216102e-08 kg/s`; this is a **physical flow imbalance** at the final state and should not be conflated with the maximum numerical conservation defect even though the magnitudes happen to be similar in this case.

The C++ `-O3 -march=native -flto -g` profiling build also passed all five CTest checks. Its Callgrind total differed from the symbol-free LTO build by only 39 instructions (`742,413,224` versus `742,413,263`), supporting comparability of the profiled code path. This does not itself prove identical native runtime performance.

## 5. Binary inspection and inlining

Symbol inspection (`nm -S --size-sort -C`) and disassembly (`objdump -d -C`) showed that separate symbols for banded solve, Jacobian assembly, Newton solve, and the friction/Reynolds helpers disappeared in both LTO builds. Fortran retained an optimized `assemble_residual_.constprop.0.isra.0` and `pipe_valid_.constprop.0`; C++ retained `pipe_sim::assemble_residual` and `pipe_sim::linepack_cz`.

The Fortran residual symbol suffixes suggest GCC interprocedural constant propagation (`constprop`) and scalar replacement of aggregates (`isra`). Their presence does **not** prove that Fortran performs less work at every corresponding source location; a missing or different symbol name is not sufficient to infer an optimization's runtime impact.

Before LTO, the Fortran banded accessor functions (`band_get_`, `band_set_`) existed as symbols, but disassembly of the optimized binary found no direct calls to them in the examined code. Therefore, **un-inlined Fortran band-accessor call overhead is not an established explanation** for the performance difference.

The C++ executable contained numerous static call sites associated with allocation, deallocation, and exception handling. Static call-site counts are not dynamic execution counts: exceptional paths may never execute. Likewise, the two static call sites for `linepack_cz` do not indicate how many times the routine runs.

## 6. Callgrind experiment

Hardware performance counters were unavailable (`kernel.perf_event_paranoid=4`), so Callgrind was used for deterministic instruction-count profiling without changing host security settings. The Callgrind runs are **not native timing benchmarks**: execution slowed to roughly 4–5 seconds under instrumentation.

| Event | Fortran LTO | C++ LTO | C++ − Fortran |
|---|---:|---:|---:|
| Total `Ir` | 660,855,731 | 742,413,263 | **81,557,532** |
| Ratio | 1.0000 | 1.1234 | +12.34% |

Callgrind's `Ir` is an instruction-reference count, not a CPU-cycle count. Different instruction mixes and memory behavior can produce different runtime ratios.

### 6.1 Math-library work matches

Exclusive instruction counts in the same major math-library functions were identical:

| Function | Fortran `Ir` | C++ `Ir` |
|---|---:|---:|
| `__ieee754_pow_fma` | 189,555,300 | 189,555,300 |
| `__ieee754_log_fma` | 78,502,700 | 78,502,700 |
| `__log10_finite` | 76,588,000 | 76,588,000 |
| `pow` wrapper | 40,208,700 | 40,208,700 |
| `log10` wrapper | 9,573,500 | 9,573,500 |

The profile recorded an edge with `calls=1914700` near `pow()`, consistent with the earlier instrumented Fortran friction-evaluation count. Identical instruction counts within these functions strongly suggest equivalent dynamic math-library work for this benchmark. They do **not** prove identical arguments or identical counts for every mathematical operation elsewhere.

This rules out **extra execution inside these measured libm routines** as the explanation for the 81.6 million excess C++ instructions. It does not establish that physics-related surrounding arithmetic, validation, or data access is equally efficient.

### 6.2 Inclusive attribution

With `callgrind_annotate --auto=no --inclusive=yes --threshold=100`, the C++ debug-symbol profile attributed approximately **530.16 million inclusive instructions (71.41%)** to residual assembly and **566.92 million (76.36%)** to source code associated with the inlined Newton path. Fortran attributed **514.73 million (77.89%)** to its optimized residual function.

**Do not sum inclusive counts**: they overlap through call relationships. Under LTO, source-line locations from inlined code are also reported under containing machine-code functions.

### 6.3 Exclusive attribution

The Fortran LTO profile's major exclusive categories included:

| Fortran location/function | Exclusive `Ir` |
|---|---:|
| `assemble_residual_.constprop.0.isra.0` | 123,817,467 |
| `MAIN__` | 102,009,736 |
| `pipe_valid_.constprop.0` | 10,144,634 |

The C++ LTO + debug-symbol profile included:

| C++ source location | Exclusive `Ir` |
|---|---:|
| `physics.cpp` attributed to residual | 80,599,443 |
| `banded_solver.cpp` attributed to `main` | 74,213,608 |
| `residual.cpp` attributed to residual | 46,765,232 |
| `newton.cpp` attributed to `main` | 35,185,985 |
| `<bits/stl_algobase.h>` attributed to `main` | 17,617,940 |
| `<bits/std_abs.h>` attributed to `main` | 13,358,079 |
| `transient_step.cpp` attributed to `main` | 8,111,604 |
| `physics.cpp` attributed to `main` | 7,770,992 |
| `integrate.cpp` attributed to `main` | 7,066,117 |
| `jacobian.cpp` attributed to `main` | 4,741,480 |

These locations are **not one-to-one counterparts** of Fortran's `MAIN__` and residual symbols. Their counts cannot directly identify the 81.6 million-instruction excess without deeper matching of the generated loops or independent kernel-level experiments.

### 6.4 Annotation/reporting lesson

Early `callgrind_annotate` commands using thresholds of `0`, `0.01`, or `1` displayed only the `__ieee754_pow_fma` function or an empty function list. `--threshold=100` exposed the complete inclusive/exclusive function summary. The profiling build contained `.debug_info`, `.debug_line`, and `.symtab` sections, allowing source-location attribution for inlined C++ code. This was a reporting issue, not evidence that Callgrind failed to collect instructions.

## 7. Source-level observations from the C++ implementation

The following observations are grounded in the current `banded_solver.cpp`, `newton.cpp`, `integrate.cpp`, and workspace declarations.

**Banded solve:** `solve_banded` uses local `get(i,j)` and `set(i,j,value)` lambdas that compute `row = kv + i - j`, check whether the band row is in range, and index `ab[j*ldab + row]`. Gaussian elimination repeatedly accesses the same pivot, uses partial pivoting, and performs back substitution. The indexing/checks may be optimized by GCC; their *dynamic* cost relative to Fortran is not yet known. Changing loop order or arithmetic grouping could alter floating-point results and requires parity testing.

**Newton solve:** `NewtonWorkspace` owns vectors for residual, trial residual, Jacobian, RHS, correction, trial state, and work state. These are sized in its constructor, and the workspace is passed by reference through the timestep/integration path. The code does **not** visibly allocate a fresh Newton workspace per iteration. It performs validation of dimensions and physical parameters, state copies, residual evaluations, Jacobian assembly, a banded solve, and a line search. Some validation repeats at different call levels; whether removing or hoisting checks is safe depends on API contracts and parity tests.

**Integration:** `IntegrationWorkspace` owns and reuses the Newton workspace and state arrays. Each step calls `transient_step`, unpacks the interleaved state, records history and linepack, evaluates the numerical conservation defect, and copies new pressure/flow back to the live state. The implementation deliberately preserves Fortran arithmetic grouping in the conservation calculation. These copies and diagnostic computations are measurable work, but are not yet demonstrated to be excessive relative to Fortran.

**No identified optimization yet:** The source offers hypotheses (band indexing, validation, copying, loop structure, compiler specialization), not a verified cause. The historical matching math-library instruction counts argued against changing the friction formula. Subsequent M13.3 optimization instead cached mathematically invariant old-time friction contributions across repeated Newton residual evaluations, reducing redundant calls without changing the physical formula.

## 8. Interpretation and limits

**Established:**

- Numerical outputs agree for the validated 8,000-step LTO benchmark, with 1,049 Newton iterations each.
- Without LTO, C++ was faster than rebuilt Fortran in the measured workload.
- Enabling LTO improved both implementations and reversed their performance ranking.
- Matched-LTO C++ is currently about 8% slower by external median process CPU time.
- C++ executes about 12.3% more Callgrind instructions overall.
- The major profiled `pow`/logarithm library functions execute identical instruction counts.

**Not established:**

- That one language is intrinsically faster.
- That a specific C++ validation check, `std::span`, `std::vector`, `std::copy`, or band-index expression causes the gap.
- That Fortran's `constprop`/`isra` symbols by themselves explain its advantage.
- That instruction-count ratios should equal elapsed-time ratios.
- That identical final printed outputs guarantee bitwise equality of all internal states.

## 9. Historical proposed experiments (superseded by M13.2b–M13.4)

1. **Freeze the baseline:** retain the current binaries or reproducible build commands, CTest results, numerical output, and benchmark summaries. Do not overwrite the validated LTO builds while experimenting.
2. **Inspect matching kernels:** compare generated instructions for the C++ and Fortran banded elimination, residual evaluation, and Newton loops. Inlining makes whole-program source mapping necessary.
3. **Use controlled microbenchmarks selectively:** benchmark matching residual, Jacobian, and banded-solve operations on identical frozen inputs; verify the isolated build reproduces relevant whole-program optimization behavior before drawing conclusions.
4. **Change one thing at a time:** begin with provably redundant bookkeeping or safe index simplifications; avoid changing physical formulas or widening numerical tolerances.
5. **Run numerical gates after each change:** five CTest checks, component parity tests, 8,000-step Fortran/C++ comparison, Newton iteration count, and linepack/mass-conservation diagnostics.
6. **Repeat controlled performance measurement:** alternating runs, external process CPU time, medians, and ideally uncertainty estimates across independent batches.
7. **Accept the original M13 C++ parity target only when** `T_CPP / T_Fortran <= 1.05` on a documented matched-build benchmark with numerical parity preserved. **This criterion was met by M13.4 against the pre-existing Fortran LTO reference; symmetric Fortran optimization remains open.**

## 10. Reproduction commands

From the repository root:

```bash
# Rebuild Fortran LTO benchmark (overwrites build/m10_cpu_banded).
make -B FC=gfortran FFLAGS_OPT="-O3 -march=native -flto" build/m10_cpu_banded

# Build native C++20 LTO executable in an isolated directory.
cmake -S . -B build/m13-cpp-lto \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -flto" \
    -DPIPE_SIM_BUILD_CPP=ON
cmake --build build/m13-cpp-lto --target pipe_sim_cpp --parallel 4
ctest --test-dir build/m13-cpp-lto --output-on-failure

# Compare final numerical outputs.
./build/m10_cpu_banded 100 15 120000
./build/m13-cpp-lto/pipe_sim_cpp 100 8000

# Callgrind; requires valgrind and callgrind_annotate.
mkdir -p benchmarks/results-m13
valgrind --tool=callgrind \
    --callgrind-out-file=benchmarks/results-m13/fortran-lto.callgrind \
    ./build/m10_cpu_banded 100 15 120000
valgrind --tool=callgrind \
    --callgrind-out-file=benchmarks/results-m13/cpp-lto.callgrind \
    ./build/m13-cpp-lto/pipe_sim_cpp 100 8000

# The threshold of 100 exposes the complete function summary.
callgrind_annotate --auto=no --inclusive=no --threshold=100 \
    benchmarks/results-m13/fortran-lto.callgrind | head -80
callgrind_annotate --auto=no --inclusive=no --threshold=100 \
    benchmarks/results-m13/cpp-lto.callgrind | head -80
```

For source attribution, use a separate C++ build with `-O3 -march=native -flto -g`, then run Callgrind on that executable. Keep the instrumented timings separate from native performance timings.

## 11. M13 optimization results (2026-10-09)

The original measurements in Sections 2–8 are retained as historical baselines, not the final status of the investigation. The C++ solver was optimized incrementally without altering the governing equations or intentionally relaxing numerical tolerances.

| Stage | Change | Evidence |
|---|---|---|
| M13.2b | Specialized direct-indexed fixed-width banded solver | Committed as `6c2730b`; dedicated solver parity (3,942 comparisons, zero error); five CTest checks passed. |
| M13.3 | Precompute old-time friction source once per Newton solve; reuse it for initial and trial residuals | Existing Fortran/C++ residual parity: 1,160 comparisons, zero relative error; additional 1,160 exact cached/uncached C++ comparisons passed; five CTest checks passed. |
| M13.4 | Add persistent `NewtonWorkspace::initial_guess`, removing per-step construction of `std::vector<double> guess(nu)` | Five CTest checks passed; full-run reference outputs unchanged; allocation removed from timestep path. |

M13.3 retains the public uncached `assemble_residual` API and adds an internal cached evaluation path. The cache is computed using the same physical-model routines and arithmetic ordering as the original old-time friction expression. The dedicated equality test checks the cases exercised by the harness; it does not establish universal bitwise identity for all possible states and compiler settings.

### 11.1 Callgrind results

| Profile | Instructions (`Ir`) |
|---|---:|
| Historical Fortran LTO | 660,855,731 |
| Original C++ LTO + symbols | 742,413,224 |
| M13.2b C++ | 709,028,666 |
| M13.3 C++ | 694,167,510 |

M13.3 eliminated **14,861,156 instructions** relative to M13.2b (**2.096%**). These are instrumented instruction counts, not native CPU time or cycle counts. No comparable M13.4 Callgrind result is yet recorded here.

### 11.2 Controlled native CPU-time benchmarks

The following figures were obtained from separate benchmark sessions; **do not subtract medians across sessions to attribute an individual optimization**. Measurements use child-process user-plus-system CPU time (`os.fork`, `os.execv`, `os.wait4`), warm-up, and 40 rounds with rotating execution order. Each session's ratio is computed from that session's medians.

| Session | Fortran LTO (ms) | Original C++ (ms) | Optimized C++ (ms) | Optimized C++ / Fortran |
|---|---:|---:|---:|---:|
| M13.2b | 132.168 | 142.909 | 139.984 | 1.05914 |
| M13.3 | 132.207 | 145.513 | 141.137 | 1.067553 |
| **M13.4** | **131.974** | **142.957** | **137.375** | **1.040921** |

In the final M13.4 session, the optimized C++ executable was **3.90% faster than the original C++ executable** (`1 - 137.375 / 142.957`) and **4.09% slower than Fortran**. The threshold was `1.05 * 131.974 = 138.573 ms`, leaving **1.198 ms** of margin. The final benchmark thus **PASSed** the stated 1.05 target for this workload and reference binary. The intermediate M13.2b and M13.3 medians should not be treated as a monotonic timing progression; session-to-session variation is evident.

No formal confidence intervals, frequency locking, or independent cross-session reproducibility analysis have yet been reported. The result establishes the recorded acceptance-test outcome, not a statistical guarantee of a fixed percentage across all workloads or systems.

### 11.3 Numerical outputs retained

For the 100-cell, 8,000-step case, the M13.4 executable produced:

| Quantity | M13.4 result |
|---|---:|
| Final outlet pressure | 7,843,657.8677442567 Pa |
| Final inlet flow | 104.99999993783898 kg/s |
| Final outlet flow | 105 kg/s |
| Final linepack | 4,798,420.4378263094 kg |
| Maximum numerical mass-conservation defect | 6.2161021219253598e-08 kg/s |
| Newton iterations | 1,049 |

These match the recorded reference outputs. The maximum numerical conservation defect must not be confused with the physical inlet/outlet flow imbalance.

## 12. Interpretation: implementation overhead versus shared mathematics

**M13.4 — allocation reuse.** The original C++ timestep constructed a new `std::vector<double>` on each call. Reusing a workspace vector removes that repeated allocation. The existing Fortran timestep declares a fixed-size `DOUBLE PRECISION GUESS(202)` local array, so it likely does not incur comparable repeated heap-allocation overhead. This is an implementation/memory-management improvement, not an algorithmic advantage over Fortran. Actual Fortran storage placement remains compiler-dependent.

**M13.3 — invariant old-time friction.** In a theta-method momentum residual, the old-time friction source depends on old pressure and old mass flow, not on the current Newton iterate. Recomputing it during each residual evaluation is mathematically redundant. The same optimization can be implemented in Fortran77. A further candidate is to precompute the entire old-time momentum contribution, including the pressure gradient, subject to numerical-parity checks and preservation of arithmetic grouping.

**M13.2b — banded indexing.** Fixed-width direct indexing reduces general indexing work in C++; analogous simplifications may or may not benefit Fortran depending on its existing band storage, compiler optimization, and generated loops. This requires inspection and measurement rather than assuming Fortran already performs the same optimization.

The observed parity therefore does **not** establish that optimized C++ is within 5% of the best attainable Fortran CPU implementation. It establishes parity with the **currently benchmarked Fortran LTO implementation**.

## 13. Agreed next phase: symmetric Fortran optimization before CUDA

Preserve the original Fortran LTO executable/source revision as a historical baseline and the validated M13.4 C++ executable/source revision as the C++ baseline. Do not silently replace either baseline. Then:

1. Inspect the Fortran residual and Newton call graph to confirm repeated old-time friction evaluation and identify suitable storage for per-timestep cached terms.
2. Implement old-time friction caching in Fortran as an isolated, reversible change. Preserve numerical evaluation order where possible, and validate residual, Newton, transient, and mass-conservation parity.
3. Consider caching the complete old-time momentum term **symmetrically in Fortran and C++**, in separate changes, and revalidate numerical equivalence.
4. Compare the Fortran and C++ banded solvers and generated code; transfer any demonstrably beneficial indexing optimizations to Fortran.
5. Benchmark original Fortran, optimized Fortran, and optimized C++ with matched compiler optimization settings, the same workload, rotating-order process-CPU measurements, and variability estimates. Record provenance and preserve all three binaries.
6. Reassess the C++/Fortran ratio against the **optimized Fortran** reference. If it exceeds 1.05, document the new gap and investigate it before describing parity against the reoptimized reference.
7. Freeze both optimized CPU baselines and their regression evidence before M14 CUDA implementation; use the best validated CPU baseline for GPU speedup claims.

The Fortran optimization phase is **planned, not completed**. No speedup or revised parity ratio against reoptimized Fortran is claimed.

## 14. Reproduction notes for the final M13.4 benchmark

The benchmarked executable was preserved as `build/m13-snapshots/pipe_sim_cpp_m13_4`. The Fortran reference was `build/m10_cpu_banded`; the original C++ comparison was `build/m13-baseline-source/pipe_sim_cpp`. Commands:

```bash
./build/m10_cpu_banded 100 15 120000
./build/m13-baseline-source/pipe_sim_cpp 100 8000
./build/m13-snapshots/pipe_sim_cpp_m13_4 100 8000
ctest --test-dir build/m13-cpp-lto --output-on-failure
```

For rigorous replication, use the original 40-round external process-CPU-time harness, not the application's internal `elapsed_s` field. The binaries and build directories are local artifacts and must be archived separately if long-term reproducibility is required. The M13.3 and M13.4 source changes were uncommitted at the time of the final reported measurements; this report does not invent a commit identifier for them.

## Conclusion

The investigation demonstrated that build configuration and source-level implementation decisions materially affect CPU performance comparisons. With matched LTO builds, the original C++20 port was slower than Fortran77. Targeted C++ banded indexing, caching of invariant old-time friction, and reuse of initial-guess storage produced a validated C++ implementation with a **1.040921×** CPU-time ratio against the existing Fortran LTO baseline for the reference workload, satisfying the original **1.05×** acceptance criterion. The original numerical outputs and all five integrated regression checks were preserved.

The next scientifically stronger comparison is **optimized Fortran versus optimized C++**, applying shared mathematical optimizations to both implementations before assessing language/compiler overhead or GPU speedup. No general conclusion about intrinsic C++ versus Fortran performance is warranted.
