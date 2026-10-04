```@meta
CurrentModule = TestParticle
```

# TestParticle.jl

TestParticle.jl is a flexible tool for tracing charged particles in electromagnetic or body force fields. It supports three-dimensional particle tracking in both relativistic and non-relativistic regimes.

The package handles field definitions in two ways:

- *Analytical Fields*: User-defined functions that calculate field values at specific spatial coordinates.

- *Numerical Fields*: Discretized fields constructed directly with coordinates or using [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) and interpolated via [FastInterpolations.jl](https://projecttorreypines.github.io/FastInterpolations.jl/stable/).

The core trajectory integration is powered by and tight to the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) ecosystem, solving the Ordinary Differential Equations (ODEs) of motion.

To accommodate different performance needs, the API provides:

- *In-place versions*: Functions ending in `!`.

- *Out-of-place versions*: Functions optimized with StaticArrays. Note that this requires the initial conditions to be passed as a static vector.

For a theoretical background on the physics involved, please refer to [Single-Particle Motions](https://henry2004y.github.io/KeyNotes/contents/single.html).

## Installation

```julia
julia> ]
pkg> add TestParticle
```

## Usage

Familiarity with the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) workflow is recommended, as TestParticle.jl builds directly upon its ecosystem.

The primary role of this package is to automate the construction of the ODE system based on Newton's second law. This allows users to focus on defining the field configurations and particle initial conditions. For practical demonstrations, please refer to the examples.

In addition to standard integrators, TestParticle.jl includes a native implementation of the Boris solver. It exposes an interface similar to DifferentialEquations.jl for ease of adoption. Further details are provided in the subsequent sections. Check more in [examples](@ref).

## Performance and Scaling

Tracing many particles can be computationally intensive. TestParticle.jl is designed to scale across multiple cores and machines using Julia's built-in parallel computing capabilities.

### Serial Performance

The following table compares the performance of TestParticle's Boris solver against standard ODE solvers from DifferentialEquations.jl for a single particle simulation.

**Benchmark Configuration:**
- Hardware: Intel Ultra 7 265K
- Task: Simulating 1 particle for 0.1 second with $dt = 1$ ns ($10^8$ steps).

| Solver | Median Time | Speedup | Allocations | Memory |
| :--- | :--- | :--- | :--- | :--- |
| TestParticle Boris | 923 ms | 1.0x | 8 | 55.12 KiB |
| ODE Tsit5 (fixed) | 15773 ms | 17x slower | 64 | 3.74 KiB |
| ODE Vern9 (fixed) | 39555 ms | 43x slower | 72 | 4.35 KiB |

### Parallel Scaling

TestParticle.jl supports both multithreading and distributed computing for ensemble simulations on CPUs. The following results were obtained on a Perlmutter 1 CPU node (2x AMD EPYC 7763). The multithreading performance is measured using `EnsembleThreads()`, while the distributed performance is measured using `EnsembleDistributed()`.

![Strong Scaling](./figures/scaling_combined.png)

### GPU Performance and Multi-Backend Scaling

TestParticle.jl provides vendor-agnostic GPU acceleration via [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl). Passing a GPU backend such as `CUDA.CUDABackend()`, `AMDGPU.ROCBackend()`, `Metal.MetalBackend()`, or `oneAPI.oneAPIBackend()` to `solve(prob, Boris(), backend)` executes the Boris pusher kernels directly on hardware accelerators.

The following benchmarks were conducted on an **NVIDIA GeForce RTX 5070** against a modern multi-core CPU. Reproducible benchmark scripts are available in [benchmark/compare_boris_backends.jl](https://github.com/henry2004y/TestParticle.jl/blob/master/benchmark/compare_boris_backends.jl).

#### 1. Analytical Uniform Field Tracing (1,000 steps)

| Number of Particles $N$ | CPU Serial (FP64) | CPU Threads (FP64) | GPU (FP64) | GPU (FP32) | Speedup (FP32 vs Serial) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 100 | 1.56 ms | 1.53 ms | 1.11 ms | 0.70 ms | 2.2x |
| 1,000 | 14.55 ms | 15.41 ms | 2.15 ms | 0.70 ms | 20.8x |
| 5,000 | 72.18 ms | 70.88 ms | 3.73 ms | 1.89 ms | 38.2x |
| 10,000 | 147.43 ms | 145.47 ms | 6.76 ms | 3.67 ms | 40.2x |

> [!NOTE]
> The table above includes host Julia overhead (allocating and assembling `ODESolution` wrapper structures). At the pure GPU kernel level without host data unpacking, an ensemble of $200,000$ particles over $1,000$ steps completes in **0.58 ms** in `Float32` ($>3.4 \times 10^{11}$ particle-steps/s), compared to **30.36 ms** in `Float64` (~52x faster kernel execution due to consumer GPU FP32:FP64 ALU architecture).

#### 2. 3D Grid Field Interpolation ($32^3$ Grid, 200 steps, Float32)

For discrete electromagnetic grids, TestParticle.jl automatically converts grid interpolations into device-native [`GPUGrid3D`](https://github.com/henry2004y/TestParticle.jl/blob/master/src/utility/gpu_grid.jl) structures with hardware register trilinear evaluation. Furthermore, passing `sort_particles = true` spatially orders particles along a 3D Morton Z-order curve prior to tracing, maximizing GPU L1/L2 cache hit rates for adjacent threads in a warp.

| Number of Particles $N$ | CPU Multithreading | GPU (Unsorted) | GPU (Morton Sorted) | Speedup (GPU vs CPU) |
| :--- | :--- | :--- | :--- | :--- |
| 500 | 2192.54 ms | 3.51 ms | 2.19 ms | ~1000x |
| 2,000 | 830.26 ms | 11.34 ms | 10.65 ms | 78.0x |
| 10,000 | 3325.21 ms | 69.23 ms | 63.11 ms | 52.7x |

### Precision Considerations: Float64 vs Float32

TestParticle.jl supports both double precision (`Float64`) and single precision (`Float32`) tracing. Supplying `Float32` coordinates, velocities, or time span automatically promotes the simulation parameters, physical constants, and equations to single precision.

#### Advantages of Float32
- **Massive GPU Speedups**: Consumer-grade GPUs (e.g., NVIDIA GeForce RTX 30/40/50 series) feature a 1:64 FP64 ALU throughput penalty. Switching to `Float32` unlocks the full FP32 compute power of the hardware, yielding up to **50x** faster pure kernel tracing.
- **Lower Memory Footprint**: Halves the required VRAM and host memory for particle trajectory states, intermediate buffers, and large 3D electromagnetic grids. This reduces memory bandwidth pressure on both CPU and GPU.

#### Trade-offs and Limitations of Float32
- **Accumulated Numerical Drift**: With 24 bits of mantissa (~7 significant decimal digits), rounding errors accumulate over long integration periods ($>10^5$ gyro-orbits), potentially leading to phase errors or slight artificial energy drift in non-integrable fields.
- **Dynamic Range / Catastrophic Cancellation**: In unnormalized SI units with large coordinate baselines (e.g. planetary magnetospheres where $r \sim 10^7$ m, but gyroradii or step increments are in centimeters or meters), `Float32` suffers from catastrophic cancellation in spatial increments.
- **Recommendation**:
  - Use **`Float64`** for high-precision single-particle tracing, long-duration orbital dynamics, and planetary-scale simulations using unnormalized SI units.
  - Use **`Float32`** for large statistical ensembles ($N \ge 10^5$), kinetic/MHD test-particle distribution studies, normalized coordinate systems (e.g. lengths normalized to ion inertial length or gyroradius), and high-throughput GPU workloads.

## Presentations

For interactive presentations and educational materials, please check out these Pluto notebooks:

- [TestParticle.jl: A New Tool for An Old Problem](https://henry2004y.github.io/pluto_playground/testparticle_202401.html)
- [基于开源工具链的测试粒子模型 (Test Particle Model Based on Open Source Toolchain)](https://henry2004y.github.io/pluto_playground/testparticle_202212.html)

## Publications

- Zijin Zhang, Anton V. Artemyev, and Vassilis Angelopoulos, 2025, "Quantification of ion scattering by solar-wind current sheets: Pitch-angle diffusion rates", Phys. Rev. E, https://doi.org/10.1103/pkzv-k5t3
- Chi Zhang, Hongyang Zhou, Chuanfei Dong, Yuki Harada, Masatoshi Yamauchi, Shaosui Xu, Hans Nilsson, et al. 2024. “Source of Drift-Dispersed Electrons in Martian Crustal Magnetic Fields”, The Astrophysical Journal, https://doi.org/10.3847/1538-4357/ad64d5.

## Acknowledgement

`TestParticle.jl` is acknowledged by citing the Zenodo DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10149789.svg)](https://doi.org/10.5281/zenodo.10149789).

Nothing can be done such easily without the support of the Julia community. We appreciate all the contributions from developers around the world.
