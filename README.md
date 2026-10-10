# Isothermal Transient Pipe Simulator

Independent research implementation of a **one-dimensional isothermal transient gas-pipeline simulator**, with a fixed-form Fortran77 reference, a validated C++20 port, and experimental CUDA acceleration.

**Status (through M14.6, October 2026):** The Fortran77 and C++20 CPU implementations are validated; CUDA component and end-to-end parity tests pass; batches of 2,048 independent GPU simulations achieve approximately **1.9× the throughput of eight CPU workers** on the tested Dell Precision 7710 / Quadro M3000M. This is a *batched throughput* result, not a single-pipeline GPU latency improvement.

## Numerical model

- Single straight, constant-diameter pipe; fixed gas temperature and composition.
- Reduced isothermal momentum equation with modular physical-property closures.
- Nodal pressure/density and staggered mass-flow unknowns.
- Conservative staggered spatial discretization, including boundary half-control volumes.
- Implicit **theta-method** time integration; theta is the time-weighting parameter.
- Prescribed inlet pressure and outlet mass flow, with transient ramp and schedule forcing in supported integrators.
- Damped Newton iteration, analytic banded Jacobian, and banded linear solver.
- Recorded linepack/inventory, boundary mass flows, Newton iterations, and mass-conservation diagnostics.

The discrete global mass-balance identity relates inventory change to theta-weighted inlet-minus-outlet mass flow. **Physical mass-flow imbalance is not numerical mass-conservation defect**: the latter measures inconsistency with the discrete inventory balance.

Details: [mathematical specification](docs/model.tex), [interface contract](docs/interface-contract.md), [milestone and validation history](docs/milestones.md).

## Implementations

| Implementation | Capabilities | Validation |
|---|---|---|
| [Fortran77](src/fortran77/) | Reference residual, Jacobian, Newton, banded solver, fixed-step and adaptive integration | Steady oracles, refinement, conservation, failure contracts, CPU optimization |
| [C++20](src/cpp/) | Native physics, residual, Jacobian, banded solver, Newton, timestep, and integration | Component and integration parity with Fortran, numerical-reference regression, performance parity |
| [CUDA](src/cuda/) | FP64 kernels, device-side Newton, persistent integration, independent-pipe batching | Component and integration parity, repeated batch benchmarks |

Test suites: [Fortran](tests/fortran77/), [C++](tests/cpp/), [CUDA](tests/cuda/). The CUDA batch implementation currently assigns **one thread per independent pipe**; Newton's method within a pipe is not yet cooperative across GPU threads.

## Milestones and reports

| Milestone | Result | Report or evidence |
|---|---|---|
| **V0–M7** | Mathematical formulation; analytic steady oracle; residual/Jacobian; damped Newton; conservative transient integration; time/space refinement; forcing, damping, and wave-response characterization | [Milestone history](docs/milestones.md) · [Model](docs/model.tex) |
| **M9** | Checked Fortran interfaces, array validation, failure and partial-progress contracts, build isolation | [Interface contract](docs/interface-contract.md) · [Milestones](docs/milestones.md) |
| **M10** | Direct-banded Jacobian and Newton solve: at N=100, dt=15 s, dense **348.155 ms**, banded **8.882 ms**, **39.20×** speedup with recorded numerical parity | [M10 performance report](benchmarks/results-m10/M10_PERFORMANCE.md) · [Baseline](benchmarks/M10_BASELINE.md) |
| **M11–M12** | Optimized Fortran reference and incremental native C++20 port with component/integration parity | [Development roadmap](docs/milestones.md) · [Cross-language investigation](docs/M13_CPP_FORTRAN_PERFORMANCE.md) |
| **M13** | Symmetric Fortran/C++ CPU optimization; final 8,000-step medians **90.998 ms Fortran**, **95.394 ms C++**; C++/Fortran ratio **1.048309**, satisfying the ≤1.05 acceptance target | [Optimization story](docs/M13_PERFORMANCE_OPTIMIZATION_STORY.md) · [Detailed performance investigation](docs/M13_CPP_FORTRAN_PERFORMANCE.md) |
| **M14.1–M14.5** | CUDA FP64 residual, Jacobian, banded solve, Newton and persistent integration with CPU/GPU parity checks | [CUDA sources](src/cuda/) · [CUDA tests](tests/cuda/) |
| **M14.6** | Parallel independent-pipe CUDA execution, multicore CPU reference, repeated timing, numerical validation; **1.950×** and **1.920×** GPU/8-worker CPU throughput at 80 and 800 steps | [M14.6 report](benchmarks/results-m14/README.md) · [CPU logs](benchmarks/results-m14/cpu-controlled.txt) · [GPU logs](benchmarks/results-m14/gpu-controlled.txt) |

**Roadmap note:** Some historical sections of [docs/milestones.md](docs/milestones.md) still carry their original *planned* labels. The later linked reports establish the actual M11–M14 work and results; the README does not reinterpret historical checklists as completed without evidence.

## Measured CPU performance (M13)

Reference case: 100 spatial cells, 15-second timestep, 8,000 steps; GCC/GFortran with `-O3 -march=native -flto`.

| Implementation | Median process CPU time |
|---|---:|
| Optimized Fortran77 | **90.998 ms** |
| Optimized C++20 | **95.394 ms** |

Both completed 1,049 Newton iterations and retained the validated numerical results. The measured ratio, **1.048309**, is workload-, hardware-, and compiler-specific, not a general ranking of programming languages. See the [M13 optimization report](docs/M13_PERFORMANCE_OPTIMIZATION_STORY.md).

## Measured GPU throughput (M14.6)

Hardware: Intel Core i7-6920HQ (4 physical / 8 logical cores) and NVIDIA Quadro M3000M (Maxwell, compute capability 5.2). Configuration: 100 cells per pipe, theta=0.65, 15-second timestep, FP64, 2,048 independent simulations; three repetitions.

| Steps per pipe | CPU, 8 workers (median) | GPU execution (median) | GPU speedup |
|---:|---:|---:|---:|
| 80 | 1,746.946 ms | 895.921 ms | **1.950×** |
| 800 | 11,850.548 ms | 6,172.766 ms | **1.920×** |

GPU execution includes host-to-device transfers, kernel, and device-to-host transfers, **but excludes allocation and CUDA runtime initialization**. CPU execution excludes the separate validation pass. Transfer overhead at this batch size was approximately 0.70% (80 steps) and 1.12% (800 steps).

Across three repetitions and the tested GPU batch sizes, **61,269,696 comparisons passed with zero failures** under the established mixed absolute/relative tolerances. Maximum observed absolute differences: 9.313e-10 (80 steps) and 1.863e-9 (800 steps). The single-pipeline GPU path remains much slower than the optimized CPU; the demonstrated advantage is for **large batches of independent simulations**.

See [full methodology, limitations, and results](benchmarks/results-m14/README.md), [raw CPU measurements](benchmarks/results-m14/cpu-controlled.txt), and [raw GPU measurements](benchmarks/results-m14/gpu-controlled.txt).

## Repository layout

| Path | Purpose |
|---|---|
| [docs/model.tex](docs/model.tex) | Governing equations, discretization, residuals, validation |
| [docs/provenance.md](docs/provenance.md) · [docs/references.bib](docs/references.bib) | Independent development and public sources |
| [docs/milestones.md](docs/milestones.md) | Roadmap and historical numerical-validation milestones |
| [docs/interface-contract.md](docs/interface-contract.md) | Checked Fortran interface, failure semantics, history contracts |
| [src/fortran77/](src/fortran77/) | Fortran reference implementation |
| [src/cpp/](src/cpp/) · [include/pipe_sim/](include/pipe_sim/) | C++20 implementation and headers |
| [src/cuda/](src/cuda/) | Experimental CUDA implementation |
| [tests/](tests/) | Regression, conservation, parity, and failure tests |
| [benchmarks/](benchmarks/) | Benchmarks, reports, measurements, provenance |

## Building and testing

The repository uses Make for Fortran workflows and CMake for C++ workflows. CUDA tests additionally require a compatible CUDA toolkit and NVIDIA GPU. Available targets and configuration depend on the toolchain.

```bash
make docs
cmake -S . -B build/cpp -DCMAKE_BUILD_TYPE=Release
cmake --build build/cpp
ctest --test-dir build/cpp --output-on-failure
```

Consult [CMakeLists.txt](CMakeLists.txt), [Makefile](Makefile), and [benchmark sources](benchmarks/) for target-specific build options. `make docs` generates `docs/model.pdf`.

## Next research directions

- **M14.7 (proposed):** profile residual evaluation, Jacobian assembly, banded linear solution, and integration on GPU; investigate cooperative Newton execution *within* one pipe, preserving numerical equations and parity.
- **CPU performance ceiling (proposed):** measure remaining initial-residual, physical-property, and timestep overheads. The [M15 CPU investigation](docs/M15_CPU_PERFORMANCE_CEILING.md) contains research hypotheses, **not achieved performance claims**.
- Future numerical assurance, physical-model/EOS extensions, adaptive surrogate studies, sensitivities, and network models remain distinct research directions.

## Independent-development principle

The implementation is derived from public governing equations and numerical literature. Prior proprietary experience may motivate questions but is **not** a source of equations, algorithms, code, data, heuristics, or tests. Numerical discrepancies are investigated causally rather than concealed by widening tolerances. Performance results are tied to their measured workloads and retained evidence.

See [provenance](docs/provenance.md) and [references](docs/references.bib).
