# M13 — C++20 / Fortran77 performance parity investigation

**Project:** Isothermal transient single-pipe simulator (`isothermal-pipe-sim`)  
**Status:** Investigation in progress; numerical parity established; matched-LTO performance criterion not yet met  
**Platform:** Dell Precision 7710, GNU GCC/GFortran 15.2, Debian Linux  
**Performance acceptance criterion:** `median(C++ CPU time) / median(Fortran CPU time) <= 1.05`  
**Reference workload:** 100 spatial cells, 15 s timestep, 8,000 steps (120,000 s simulated)

## Executive summary

The project ported a validated fixed-form Fortran77 transient gas-pipeline solver to native C++20, retaining the numerical formulation and testing each translated component against its Fortran counterpart. Initial measurements suggested that C++ was faster than Fortran. Controlled recompilation revealed that this result depended strongly on **link-time optimization (LTO)**: C++ led in matched non-LTO builds, whereas Fortran led when both implementations were compiled with `-O3 -march=native -flto`.

The latest alternating-run, externally measured process-CPU-time experiment found medians of **131.778 ms (Fortran LTO)** and **142.291 ms (C++ LTO)**, a C++/Fortran ratio of **1.07978**. The M13 acceptance ratio is 1.05; at that measured Fortran baseline, C++ must reach approximately **138.367 ms**, a reduction of **3.924 ms** (about 2.76% of current C++ time). The ratio is workload- and machine-specific; this is not a general ranking of languages.

Crucially, both LTO executables produced matching reported final outlet pressure, linepack, Newton iteration count, and maximum numerical mass-conservation defect. Callgrind reported **660,855,731 instructions** for Fortran and **742,413,263** for C++, a difference of **81,557,532** (12.34% relative to Fortran). The dominant measured libm routines had **identical exclusive instruction counts**, so the excess is not explained by additional execution *inside those routines*. Attribution of the remaining difference to a specific C++ kernel or compiler transformation is ongoing.

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

**No identified optimization yet:** The source offers hypotheses (band indexing, validation, copying, loop structure, compiler specialization), not a verified cause. The matching math-library instruction counts make wholesale changes to the friction formula an especially poor first experiment.

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

## 9. Proposed next experiments

1. **Freeze the baseline:** retain the current binaries or reproducible build commands, CTest results, numerical output, and benchmark summaries. Do not overwrite the validated LTO builds while experimenting.
2. **Inspect matching kernels:** compare generated instructions for the C++ and Fortran banded elimination, residual evaluation, and Newton loops. Inlining makes whole-program source mapping necessary.
3. **Use controlled microbenchmarks selectively:** benchmark matching residual, Jacobian, and banded-solve operations on identical frozen inputs; verify the isolated build reproduces relevant whole-program optimization behavior before drawing conclusions.
4. **Change one thing at a time:** begin with provably redundant bookkeeping or safe index simplifications; avoid changing physical formulas or widening numerical tolerances.
5. **Run numerical gates after each change:** five CTest checks, component parity tests, 8,000-step Fortran/C++ comparison, Newton iteration count, and linepack/mass-conservation diagnostics.
6. **Repeat controlled performance measurement:** alternating runs, external process CPU time, medians, and ideally uncertainty estimates across independent batches.
7. **Accept M13 only when** `T_CPP / T_Fortran <= 1.05` on a documented matched-build benchmark with numerical parity preserved.

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

## Conclusion

This investigation illustrates why scientific-software performance claims must control **algorithm, compiler, optimization mode, workload, numerical results, and measurement method**. The C++20 port initially appeared faster than Fortran77. Once both implementations used LTO, Fortran became faster, even though their expensive mathematical library work and final numerical outputs matched. The remaining gap is modest but measurable; identifying its source is the outstanding M13 engineering task. The evidence supports targeted investigation of generated numerical loops and bookkeeping, **not** a general conclusion about C++ versus Fortran performance.
