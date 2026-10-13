========================================================================================
Boris Pusher Benchmark: CPU Serial vs CPU Threads vs GPU Backends
CPU Threads: 128
GPU Backend: CUDA
GPU Device:  NVIDIA A100-SXM4-40GB
Timing: median of 3 runs after one warmup; save_everystep = false
========================================================================================

[1] Analytical Uniform Field Tracing (1000 steps)
----------------------------------------------------------------------------------------
Particles | CPU Ser    | CPU Thr    | KA CPU     | GPU (FP64) | GPU (FP32)
----------------------------------------------------------------------------------------
1000      |   35.86 ms |    1.88 ms |    1.33 ms |    1.48 ms |    1.32 ms
10000     |  337.73 ms |   10.05 ms |    3.40 ms |    2.06 ms |    1.63 ms
100000    | 3642.52 ms |  262.42 ms |   24.39 ms |   10.02 ms |    4.93 ms

  Speedup relative to CPU Serial:
Particles | CPU Thr       | KA CPU        | GPU (FP64)    | GPU (FP32)   
--------------------------------------------------------------------------
1000      | 19.08x        | 26.92x        | 24.30x        | 27.15x       
10000     | 33.62x        | 99.45x        | 164.16x       | 207.15x      
100000    | 13.88x        | 149.34x       | 363.63x       | 738.94x      

  Raw output (bulk `u0` matrix, no ODESolution per particle):
Particles | KA CPU        | GPU (FP64)    | GPU (FP32)   
----------------------------------------------------------------------------------------
1000      |       0.78 ms |       1.51 ms |       0.57 ms
10000     |       2.30 ms |       0.71 ms |       0.57 ms
100000    |      20.23 ms |       2.33 ms |       1.22 ms

  Speedup relative to KA CPU, raw output:
Particles | GPU (FP64)    | GPU (FP32)   
------------------------------------------
1000      | 0.52x         | 1.37x        
10000     | 3.25x         | 4.04x        
100000    | 8.67x         | 16.63x       

  Effective throughput at 100000 particles (1000 steps):
    CPU serial     27.5 | CPU threads    381.1 | KA CPU   4099.8 | GPU FP64   9983.0 | GPU FP32  20286.5 M particle-steps/s

[2] 3D Grid Field Interpolation (32³ Grid, 200 steps, Float32)
----------------------------------------------------------------------------------------
Particles | CPU Thr         | KA CPU          | GPU            
----------------------------------------------------------------------------------------
2000      |         2.60 ms |         1.73 ms |         1.57 ms
10000     |         9.85 ms |         3.39 ms |         1.73 ms
50000     |        85.47 ms |        11.06 ms |         3.81 ms

  `CPU Thr` builds one ODESolution per particle, which with grid fields
  costs far more than the pusher itself; `KA CPU` is the fair CPU baseline.

  Speedup relative to CPU Threads:
Particles | KA CPU        | GPU          
------------------------------------------
2000      | 1.50x         | 1.65x        
10000     | 2.90x         | 5.70x        
50000     | 7.73x         | 22.41x       

  Raw output (bulk `u0` matrix, no ODESolution per particle):
Particles | KA CPU          | GPU            
----------------------------------------------------------------------------------------
2000      |         1.04 ms |         0.78 ms
10000     |         2.56 ms |         0.78 ms
50000     |        10.29 ms |         1.05 ms

  Speedup relative to KA CPU, raw output:
Particles | GPU          
--------------------------
2000      | 1.34x        
10000     | 3.27x        
50000     | 9.82x        

[3] GPU-only large-scale throughput (1000 steps)
----------------------------------------------------------------------------------------
Particles | GPU (FP64)    | GPU (FP32)    | FP32 raw      | FP32 M/s     | Raw M/s     
----------------------------------------------------------------------------------------
100000    |       9.14 ms |       4.76 ms |       1.33 ms |      20990.6 |      75463.3
500000    |     134.28 ms |      21.15 ms |       4.86 ms |      23638.6 |     102834.8
1000000   |     341.07 ms |     162.05 ms |       9.66 ms |       6171.1 |     103501.9
5000000   |    1883.00 ms |     829.77 ms |      37.43 ms |       6025.8 |     133585.0

  Speedup relative to GPU (FP64):
Particles | GPU (FP32)    | FP32 raw     
------------------------------------------
100000    | 1.92x         | 6.90x        
500000    | 6.35x         | 27.62x       
1000000   | 2.10x         | 35.30x       
5000000   | 2.27x         | 50.31x       
  `FP32 raw` uses a bulk `u0` matrix and `raw_output = true`, so no
  ODESolution is built. It changes the host bookkeeping only.

[4] GPU kernel isolation, analytical field (1000000 particles, step sweep)
----------------------------------------------------------------------------------------
Steps     | FP64          | FP32          | FP64 raw      | FP32 raw     
----------------------------------------------------------------------------------------
100       |     328.21 ms |     155.14 ms |       9.23 ms |       5.92 ms
1000      |     336.31 ms |     160.46 ms |      16.03 ms |       8.91 ms
10000     |     397.34 ms |     187.91 ms |      71.37 ms |      39.87 ms
50000     |     592.28 ms |     283.87 ms |     270.75 ms |     139.08 ms
200000    |    1377.83 ms |     677.32 ms |    1059.02 ms |     530.25 ms

  Speedup relative to analytical FP64:
Steps     | FP32          | FP32 raw     
------------------------------------------
100       | 2.12x         | 55.47x       
1000      | 2.10x         | 37.73x       
10000     | 2.11x         | 9.97x        
50000     | 2.09x         | 4.26x        
200000    | 2.03x         | 2.60x        

  Fit of t = a + b*steps at 1000000 particles
    case             variant   | host cost a   | kernel per step b | kernel M steps/s
    ----------------------------------------------------------------------------
    analytical FP64  baseline  |     333.82 ms |         0.0052 ms |        191584.4
    analytical FP64  raw       |      12.00 ms |         0.0052 ms |        191088.6
    analytical FP32  baseline  |     157.32 ms |         0.0026 ms |        385058.2
    analytical FP32  raw       |       8.51 ms |         0.0026 ms |        383132.2
    `a` is the host cost (initialization and ODESolution assembly);
    `b` is the pusher kernel alone. A negative `a` means the model
    does not hold and the fit cannot be trusted. `b` must agree
    between variants, since the raw path changes the host side only.

[5] GPU kernel isolation, 32³ grid field (1000000 particles, step sweep)
----------------------------------------------------------------------------------------
Steps     | grid FP32     | grid FP32 raw
----------------------------------------------------------------------------------------
10        |      92.93 ms |       6.27 ms
100       |      95.31 ms |       9.09 ms
1000      |     122.53 ms |      34.81 ms
10000     |     276.26 ms |     192.26 ms

  Fit of t = a + b*steps at 1000000 particles
    case             variant   | host cost a   | kernel per step b | kernel M steps/s
    ----------------------------------------------------------------------------
    grid FP32        baseline  |      96.67 ms |         0.0180 ms |         55451.2
    grid FP32        raw       |       9.73 ms |         0.0183 ms |         54592.0
    `a` is the host cost (initialization and ODESolution assembly);
    `b` is the pusher kernel alone. A negative `a` means the model
    does not hold and the fit cannot be trusted. `b` must agree
    between variants, since the raw path changes the host side only.

  Grid interpolation costs 7.0x the analytical field per step

========================================================================================