========================================================================================
Boris Pusher Benchmark: CPU Serial vs CPU Threads vs GPU Backends
CPU Threads: 20
GPU Backend: CUDA
GPU Device:  NVIDIA GeForce RTX 5070
Timing: median of 3 runs after one warmup; save_everystep = false
========================================================================================

[1] Analytical Uniform Field Tracing (1000 steps)
----------------------------------------------------------------------------------------
Particles | CPU Ser    | CPU Thr    | KA CPU     | GPU (FP64) | GPU (FP32)
----------------------------------------------------------------------------------------
1000      |   14.15 ms |    1.62 ms |    0.79 ms |    2.27 ms |    0.52 ms
10000     |  167.31 ms |   12.89 ms |    5.76 ms |    3.09 ms |    0.85 ms
100000    | 1521.67 ms |  166.30 ms |   43.22 ms |   23.40 ms |    5.30 ms

  Speedup relative to CPU Serial:
Particles | CPU Thr       | KA CPU        | GPU (FP64)    | GPU (FP32)
--------------------------------------------------------------------------
1000      | 8.71x         | 17.96x        | 6.24x         | 27.42x
10000     | 12.98x        | 29.06x        | 54.21x        | 196.86x
100000    | 9.15x         | 35.21x        | 65.04x        | 286.89x

  Raw output (bulk `u0` matrix, no ODESolution per particle):
Particles | KA CPU        | GPU (FP64)    | GPU (FP32)
----------------------------------------------------------------------------------------
1000      |       0.81 ms |       1.87 ms |       0.45 ms
10000     |       5.26 ms |       2.01 ms |       0.60 ms
100000    |      37.69 ms |      13.29 ms |       1.22 ms

  Speedup relative to KA CPU, raw output:
Particles | GPU (FP64)    | GPU (FP32)
------------------------------------------
1000      | 0.43x         | 1.78x
10000     | 2.62x         | 8.83x
100000    | 2.84x         | 30.82x

  Effective throughput at 100000 particles (1000 steps):
    CPU serial     65.7 | CPU threads    601.3 | KA CPU   2313.9 | GPU FP64   4274.0 | GPU FP32  18853.7 M particle-steps/s

[2] 3D Grid Field Interpolation (32³ Grid, 200 steps, Float32)
----------------------------------------------------------------------------------------
Particles | CPU Thr         | KA CPU          | GPU
----------------------------------------------------------------------------------------
2000      |         3.39 ms |         1.67 ms |         0.95 ms
10000     |        11.29 ms |         6.34 ms |         1.64 ms
50000     |        79.71 ms |        23.99 ms |         4.40 ms

  `CPU Thr` builds one ODESolution per particle, which with grid fields
  costs far more than the pusher itself; `KA CPU` is the fair CPU baseline.

  Speedup relative to CPU Threads:
Particles | KA CPU        | GPU
------------------------------------------
2000      | 2.04x         | 3.57x
10000     | 1.78x         | 6.90x
50000     | 3.32x         | 18.12x

  Raw output (bulk `u0` matrix, no ODESolution per particle):
Particles | KA CPU          | GPU
----------------------------------------------------------------------------------------
2000      |         1.42 ms |         0.75 ms
10000     |         5.79 ms |         0.74 ms
50000     |        22.85 ms |         1.26 ms

  Speedup relative to KA CPU, raw output:
Particles | GPU
--------------------------
2000      | 1.89x
10000     | 7.84x
50000     | 18.13x

[3] GPU-only large-scale throughput (1000 steps)
----------------------------------------------------------------------------------------
Particles | GPU (FP64)    | GPU (FP32)    | FP32 raw      | FP32 M/s     | Raw M/s
----------------------------------------------------------------------------------------
100000    |      22.42 ms |       5.81 ms |       1.23 ms |      17220.6 |      81453.1
500000    |     144.52 ms |      47.72 ms |       3.92 ms |      10478.3 |     127674.8
1000000   |     292.05 ms |      89.99 ms |       8.62 ms |      11112.9 |     115973.0
5000000   |    1568.81 ms |     632.52 ms |      43.10 ms |       7904.8 |     116007.7

  Speedup relative to GPU (FP64):
Particles | GPU (FP32)    | FP32 raw
------------------------------------------
100000    | 3.86x         | 18.26x
500000    | 3.03x         | 36.90x
1000000   | 3.25x         | 33.87x
5000000   | 2.48x         | 36.40x
  `FP32 raw` uses a bulk `u0` matrix and `raw_output = true`, so no
  ODESolution is built. It changes the host bookkeeping only.

[4] GPU kernel isolation, analytical field (1000000 particles, step sweep)
----------------------------------------------------------------------------------------
Steps     | FP64          | FP32          | FP64 raw      | FP32 raw
----------------------------------------------------------------------------------------
100       |     199.29 ms |     100.40 ms |      24.31 ms |       7.12 ms
1000      |     283.78 ms |      98.93 ms |     119.55 ms |       8.78 ms
10000     |    1251.91 ms |     125.33 ms |    1075.63 ms |      23.85 ms
50000     |    5500.48 ms |     180.37 ms |    5318.68 ms |      90.89 ms

  Speedup relative to analytical FP64:
Steps     | FP32          | FP32 raw
------------------------------------------
100       | 1.99x         | 27.99x
1000      | 2.87x         | 32.32x
10000     | 9.99x         | 52.50x
50000     | 30.50x        | 60.51x

  Fit of t = a + b*steps at 1000000 particles
    case             variant   | host cost a   | kernel per step b | kernel M steps/s
    ----------------------------------------------------------------------------
    analytical FP64  baseline  |     184.69 ms |         0.1063 ms |          9404.8
    analytical FP64  raw       |      13.88 ms |         0.1061 ms |          9425.2
    analytical FP32  baseline  |     101.83 ms |         0.0016 ms |        625405.4
    analytical FP32  raw       |       7.04 ms |         0.0017 ms |        596209.8
    `a` is the host cost (initialization and ODESolution assembly);
    `b` is the pusher kernel alone. A negative `a` means the model
    does not hold and the fit cannot be trusted. `b` must agree
    between variants, since the raw path changes the host side only.

[5] GPU kernel isolation, 32³ grid field (1000000 particles, step sweep)
----------------------------------------------------------------------------------------
Steps     | grid FP32     | grid FP32 raw
----------------------------------------------------------------------------------------
10        |     107.27 ms |       7.02 ms
100       |     104.02 ms |       8.26 ms
1000      |     121.19 ms |      30.94 ms
10000     |     311.17 ms |     217.95 ms

  Fit of t = a + b*steps at 1000000 particles
    case             variant   | host cost a   | kernel per step b | kernel M steps/s
    ----------------------------------------------------------------------------
    grid FP32        baseline  |     103.24 ms |         0.0208 ms |         48160.5
    grid FP32        raw       |       7.54 ms |         0.0211 ms |         47476.1
    `a` is the host cost (initialization and ODESolution assembly);
    `b` is the pusher kernel alone. A negative `a` means the model
    does not hold and the fit cannot be trusted. `b` must agree
    between variants, since the raw path changes the host side only.

  Grid interpolation costs 12.6x the analytical field per step

========================================================================================
Tips:
 - To run on GPUs (e.g. CUDA, AMDGPU, Metal):
     julia --project=test -t auto benchmark/compare_boris_backends.jl CUDA
 - To run with multiple CPU threads:
     julia --project=test -t auto benchmark/compare_boris_backends.jl
 - Sections [1] and [2] report the default API first, then repeat the
   EnsembleKernel paths with a bulk `u0` matrix and `raw_output = true`.
 - Thread scaling and multi-node scaling live in run_scaling_threads.jl and
   run_scaling_multinode.jl.
========================================================================================