# M14.7 — Newton Operation Census and GPU Resource Characterization

## CPU Newton census

Workload: 8,000 timesteps.

| Measurement | Value |
|---|---:|
| Total Newton iterations | 1,049 |
| Mean iterations per timestep | 0.131125 |
| Maximum iterations | 2 |
| Zero-iteration timesteps | 7,142 |
| One-iteration timesteps | 667 |
| Two-iteration timesteps | 191 |
| Last nonzero-iteration timestep | 858 |

The CPU diagnostic passed.

## CUDA resource characterization

Target: NVIDIA Quadro M3000M, compute capability 5.2.

| Resource | Persistent integration | Newton-only |
|---|---:|---:|
| Registers per thread | 131 | 112 |
| Stack frame (bytes/thread) | 22,608 | 20,192 |
| Compiler-reported spill stores | 0 | 0 |
| Compiler-reported spill loads | 0 | 0 |

Large stack frames indicate substantial thread-local storage
requirements but do not establish runtime memory traffic.

## Nsight Systems

Workload: 80 timesteps per pipeline.

| Batch size | Blocks | Threads/block | Kernel duration (ms) |
|---|---:|---:|---:|
| 512 | 16 | 32 | 809.389 |
| 1,024 | 32 | 32 | 760.778 |
| 2,048 | 64 | 32 | 893.709 |

All three batch parity checks passed.

These timings were collected under Nsight Systems and should
not be treated as unprofiled throughput measurements.

Nsight Compute hardware-counter collection was unavailable
because of ERR_NVGPUCTRPERM.

## Next milestone

M14.7c: optional device-side timing of Newton setup,
residual evaluation, Jacobian assembly, banded solve,
line search, and integration overhead.
