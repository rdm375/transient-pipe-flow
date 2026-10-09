# M14.6 — Parallel CUDA Pipeline Simulation Benchmarks

## Hardware

CPU: Intel Core i7-6920HQ
- 4 physical cores
- 8 logical processors

GPU: NVIDIA Quadro M3000M
- Maxwell GM204
- Compute capability 5.2
- 4 GiB device memory

## Configuration

- 100 spatial cells per pipeline
- Double-precision arithmetic
- Theta-method time integration
- Theta = 0.65
- Timestep = 15 seconds
- Batch size = 2048 independent simulations
- Workloads = 80 and 800 timesteps
- Three repetitions per configuration

CPU:
- GCC 15
- C++20
- -O3 -march=native -flto
- 1, 2, 4, and 8 workers

GPU:
- CUDA 12.4
- sm_52
- -O3 --fmad=false
- One CUDA thread per independent pipeline

## Median Performance

| Steps | CPU 1 worker (ms) | CPU 4 workers (ms) | CPU 8 workers (ms) | GPU execution (ms) |
|---|---:|---:|---:|---:|
| 80 | 8637.859 | 2530.396 | 1746.946 | 895.921 |
| 800 | 58168.918 | 16986.292 | 11850.548 | 6172.766 |

GPU execution includes host-to-device transfers, kernel
execution, and device-to-host transfers. It excludes device
allocation and CUDA runtime initialization.

CPU execution excludes numerical validation.

## GPU Speedup Relative to Eight CPU Workers

80 steps:

    1746.946 / 895.921 = 1.950x

800 steps:

    11850.548 / 6172.766 = 1.920x

## GPU Transfer Overhead

At batch size 2048:

| Steps | H2D (ms) | Kernel wall time (ms) | D2H (ms) |
|---|---:|---:|---:|
| 80 | 0.494 | 889.689 | 5.800 |
| 800 | 0.438 | 6103.825 | 68.502 |

Transfers account for approximately 0.70% and 1.12%
of execution time, respectively.

## Numerical Validation

All CPU benchmark configurations passed exact
CPU-to-CPU validation.

All GPU configurations passed the established
mixed absolute/relative numerical tolerances.

Across three repetitions:

- 80 steps: 7,411,968 comparisons, zero failures
- 800 steps: 53,857,728 comparisons, zero failures

Total: 61,269,696 comparisons, zero failures.

Maximum observed absolute differences:

- 80 steps: 9.313e-10
- 800 steps: 1.863e-9

## Interpretation

The GPU achieves approximately 1.9x the throughput
of eight CPU workers for sufficiently large batches
of independent transient pipeline simulations.

The performance advantage is primarily computational;
host/device transfers contribute little to total time.

The current CUDA implementation assigns one thread
to each independent pipeline. It does not parallelize
Newton's method within an individual pipeline.

Consequently, these results demonstrate batched
throughput acceleration, not single-pipeline latency
acceleration.

## Limitations

- Results apply to the tested hardware and workloads.
- Three repetitions were performed.
- CPU and GPU measurements were collected separately.
- GPU execution excludes device allocation.
- CPU timing includes state reset and thread launch overhead.
- The precise batch-size crossover was not measured.
- No cooperative GPU Newton solver was evaluated.

## Raw Measurements

- cpu-controlled.txt
- gpu-controlled.txt
