# Features Walkthrough

`TestParticle.jl` provides a range of features to simplify and enhance your particle tracing simulations. This section provides a walkthrough of these capabilities with links to specific demos.

## Summary of Features

The table below summarizes the key features and tracing capabilities in `TestParticle.jl`:

| Category | Capability | Description | Reference |
| :--- | :--- | :--- | :--- |
| **Solvers & Dynamics** | [Boris Pusher](@ref Boris-Pusher) | Symplectic energy-conserving pushers (Standard, Multistep, Hyper, Adaptive) | [Tutorial](@ref Boris-Pusher), [Demo](@ref Boris-Method) |
| | SciML Solvers | High-order adaptive and symplectic ODE solvers from `OrdinaryDiffEq.jl` | [Numerical Methods](@ref) |
| | Guiding Center (GC) | Reduced-order drift dynamics averaging over gyromotion | [Guiding Center](@ref), [Guiding Center Approximation](@ref) |
| | Hybrid Solver | Adaptive runtime switching between full orbit and guiding center | [Hybrid GC-FO Solver](@ref) |
| | Canonical Form | Relativistic/non-relativistic canonical $(x, p)$ phase space with auto-diff | [Canonical Tracing Form](@ref) |
| | Fieldline Tracing | Traces magnetic field lines in 3D fields | [Tracing Magnetic Field Lines](@ref) |
| | Relativistic Dynamics | Relativistic equations of motion for high-energy particles | [Relativistic vs Non-relativistic Tracing](@ref) |
| | Backward Tracing | Reverse-time integration for source identification | [Backwards in Time](#tracing-backwards-in-time) |
| **Callbacks & Boundary** | `isoutside` Predicate | Domain boundary check with step rollback (CPU & GPU) | [Boundary Checks](#callbacks-and-boundary-checks) |
| | `TerminateOutside` | SciML discrete callback wrapper for boundary termination | [Boundary Checks](#callbacks-and-boundary-checks) |
| | SciML Callbacks | Continuous/discrete event handling and step size limiters | [Numerical Methods](@ref) |
| **Parallel & Hardware** | Multithreading | CPU shared-memory parallelism (`EnsembleThreads`) | [Multithreaded Parallelization](@ref) |
| | Distributed & Hybrid | Multi-process and multi-node clusters (`EnsembleDistributed`, `EnsembleSplitThreads`) | [Distributed Parallelization](@ref) |
| | Native GPU Kernels | Backend-agnostic GPU execution (CUDA, AMDGPU, oneAPI, Metal) | [GPU Ensemble Tracing](@ref) |
| **Output & Diagnostics** | `save_fields` & `save_work` | Inline trajectory logging of fields ($\mathbf{E}, \mathbf{B}$) and work rates | [Extra Diagnostics](@ref) |
| | `saveat` & Dense Output | Explicit output times and linear step interpolation | [Boris Pusher](@ref Boris-Pusher) |
| | `raw_output` | High-throughput dense array output bypassing solution wrapping on GPU | [GPU Ensemble Tracing](@ref) |
| **Fields & Units** | Gridded Fields | 1D/2D/3D regular, rectilinear, and structured grid interpolation | [Field Interpolation](@ref) |
| | Spherical Device Grids | GPU-optimized trilinear spherical grid interpolator (`GPUSphericalGrid`) | [GPU Ensemble Tracing](@ref) |
| | Dimensionless Units | Normalized scaling by gyrofrequency, gyroradius, etc. | [Dimensionless Units](@ref Dimensionless-Units-and-Normalization) |
| | Unitful Support | Physical units tracking via `Unitful.jl` | [Unitful Integration](@ref Unit-support) |
| | Velocity Distributions | Particle sampling (Maxwellian, Kappa) from distribution functions | [Velocity Sampling](@ref) |

## Unit Conversions and Scaling

By default, SI units are used. However, you can define your own unit systems by setting reference scales for mass, charge, magnetic field, length, and time.

- **Check out the demos on unit conversions**:
    - [Dimensionless Units](@ref Dimensionless-Units-and-Normalization)
    - [Unitful Integration](@ref Unit-support)

## Tracing Backwards in Time

Tracing a particle backwards is as simple as setting a time span where the end time is less than the start time, e.g., `tspan = (0.0, -1.0)`. For fixed time step methods, you also need to set `dt` to be negative (e.g. standard Boris pusher). This is particularly useful for source identification in plasma simulations.

## Callbacks and Boundary Checks

In particle tracing simulations, it is frequently necessary to stop integration when a
particle leaves the simulation box, collides with a boundary or body, or encounters a
specific physical condition. `TestParticle.jl` provides two primary mechanisms for handling
boundary checks and callbacks:

### 1. The `isoutside` Predicate

Used by the **Boris family** (both host CPU and native GPU `KernelAbstractions` backends),
**Guiding Center** (`TraceGCProblem`), and **Hybrid** (`TraceHybridProblem`) solvers.

- **Signature**: `isoutside(u, p, t) -> Bool`, where `u` is the phase-space state vector,
  `p` is the parameter container, and `t` is current time.
- **Behavior**: When `isoutside` evaluates to `true`, the solver rolls back the step that
  crossed the boundary so the trajectory ends at the last state inside the domain, and marks
  the solution with `ReturnCode.Terminated`.
- **GPU Portability**: On device backends (e.g., `CUDA.CUDABackend()`), `isoutside` is
  passed as an argument to the device kernel and executed directly on the GPU. It must
  therefore be an `isbits` function that avoids dynamic allocations.
- **High-throughput ensembles**: In `raw_output = true` mode, the `terminated` boolean array
  marks which particles exited early, and unused time slots repeat the final inside state.

```julia
using TestParticle, StaticArrays
using TestParticle: Kappa, Maxwellian

# Stop once the particle exits the region |x| <= 10 m
isoutside(u, p, t) = abs(u[1]) > 10.0

# Native Boris solver on CPU
sol = solve(prob, Boris(); dt = 1.0e-9, isoutside)

# Native Boris solver on GPU
# using CUDA
# sol = solve(prob, Boris(), CUDABackend(); dt = 1.0e-9, trajectories = 100, isoutside)
```

### 2. `TerminateOutside` and SciML Callbacks

Used by standard **SciML ODE solvers** (`OrdinaryDiffEq.jl`, e.g., `Tsit5()`, `Vern9()`).

- **`TerminateOutside(isoutside)`**: Convenience wrapper that constructs a SciML
  `DiscreteCallback` triggering `terminate!(integrator)` when `isoutside(u, p, t)` holds.
- **SciML Callback Interface**: Standard `ContinuousCallback`, `DiscreteCallback`, or
  `CallbackSet` can be passed via `callback = ...` for generic event handling, root-finding,
  or enforcing maximum step size limits based on gyroperiods.

```julia
using TestParticle, OrdinaryDiffEq
using TestParticle: Kappa, Maxwellian

# Domain termination for SciML ODE solvers
isoutside(u, p, t) = hypot(u[1], u[2], u[3]) > 10.0
callback = TerminateOutside(isoutside)

sol = solve(prob, Tsit5(); callback)
```

### Comparison Matrix

| Feature | Solvers | Hardware | Behavior |
| :--- | :--- | :--- | :--- |
| **`isoutside`** | `Boris()`, `MultistepBoris`, GC, Hybrid | CPU & GPU | Step rollback to boundary; sets `ReturnCode.Terminated` |
| **`TerminateOutside`** | SciML ODE solvers (`Tsit5`, `Vern9`, etc.) | CPU | Halts integrator via `terminate!` |
| **SciML Callbacks** | SciML ODE solvers | CPU | Arbitrary event handling, root-finding, stepsize limiters |

For detailed tutorials and examples, see:
- [Boundary check](@ref) in Numerical Methods
- [Boundary Checks on Device](@ref) in Boris Pusher
- [Boundary Checks with isoutside](@ref) in GPU Ensemble Tracing

## Multiple Particles (Ensemble Simulations)

Tracing multiple particles simultaneously is efficiently handled through ensemble simulations. `TestParticle.jl` supports several parallelization strategies to leverage available hardware (following [Ensemble Algorithms](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/#EnsembleAlgorithms)):

- **EnsembleSerial**: Sequential execution on a single core.
- **EnsembleThreads**: Shared-memory parallelization using multiple threads. This is the most efficient choice for local multi-core machines.
- **EnsembleDistributed**: Distributed-memory parallelization across multiple processes or nodes. Useful for large-scale simulations on a cluster.
- **EnsembleSplitThreads**: A hybrid strategy that combines distributed processes with multi-threading on each process. This is ideal for HPC clusters where you want to use multiple nodes and all cores within each node.

### Work Distribution (`batch_size`)

When using distributed parallelization (`EnsembleDistributed` or `EnsembleSplitThreads`), you can control the granularity of task distribution using the `batch_size` keyword argument in `solve`. This parameter determines how many trajectories are sent to a worker in a single chunk.

By default, `batch_size` is automatically optimized for the current environment:
- For distributed solvers, it defaults to `max(1, trajectories ÷ nworkers())`.
- For serial and threaded solvers, it defaults to 1.

```julia
# Explicitly setting batch_size for distributed simulation
sols = solve(prob, EnsembleSplitThreads(); trajectories=10000, batch_size=1000)
```

### Efficiency Considerations

When tracing in numeric fields, memory efficiency is critical. By default, ensure that large field data is shared as a reference rather than being copied for each trajectory. 

- **Example**: [Ensemble Tracing](@ref Ensemble-Tracing)

!!! warning
    When defining a `prob_func` for an `EnsembleProblem`, avoid in-place mutation of the `ODEProblem` (e.g., `prob.u0 = ...`), especially when using `StaticArrays`. SciMLBase v2.0+ may issue a warning about "Mutation of ODEProblem detected". Instead, use `remake(prob; u0=...)` to create a new problem instance with updated parameters.

## Parallelization

Tracing many independent particles is an *embarrassingly parallel* workload: each trajectory is computed without communication with the others. `TestParticle.jl` exposes the same ensemble parallelization interface for both the SciML solvers (via `OrdinaryDiffEq`) and the native [Boris Method](@ref), so switching from a serial run to a parallel one only requires changing a single argument.

| Algorithm | Description |
| :--- | :--- |
| `EnsembleSerial()` | Trajectories are solved one after another (no parallelism). |
| `EnsembleThreads()` | Trajectories are split across the available Julia threads within a single process. |
| `EnsembleDistributed()` | Trajectories are distributed across worker processes with `pmap` (multiple machines / NUMA nodes). |
| `EnsembleSplitThreads()` | Hybrid: work is first split across worker processes, then multithreaded within each process. |

The two execution models differ enough that they are documented separately:

- [Multithreaded Parallelization](@ref) (`EnsembleSerial` / `EnsembleThreads`) runs entirely inside one Julia session — no extra setup is needed, and the achievable speedup depends on `Threads.nthreads()`.
- [Distributed Parallelization](@ref) (`EnsembleDistributed` / `EnsembleSplitThreads`) requires adding worker processes up front and loading the required packages on every worker, but scales beyond a single machine.

For GPU-based massively parallel tracing, see the [GPU Ensemble Tracing](@ref) example.

## [Boris Pusher](@id features-boris-pusher)

The native Boris pusher is optimized for performance and energy conservation in magnetic fields. `TestParticle.jl` provides a whole family of Boris solvers — standard, Multistep, Hyper, and Adaptive — described in [Boris Pusher](@ref Boris-Pusher). `MultistepBoris6(; n = 4)` is the recommended default for most work.

- **Tutorial**: [Boris Pusher](@ref Boris-Pusher)
- **Example**: [Boris Method](@ref)

## Numerical Field Interpolations

`TestParticle.jl` supports tracing in numerical fields defined on various grids (Cartesian, Rectilinear, Structured).

- **Example**: [Field Interpolation](@ref)

## Real-world Applications

Explore how `TestParticle.jl` can be applied to complex physical systems:

- [Shock Acceleration](@ref Shock)
- [Radiation Effects](@ref Radiation)
- [Cosmic Ray Tracing](@ref Cosmic-Ray)

## More Demos

For a full list of capabilities, browse the **Examples** section in the sidebar:

- **Drifts**: [ExB Drift](@ref), [Curvature and Grad-B Drifts](@ref),
  [Magnetic Drift and Energy Partition](@ref),
  [Relation between the Magnetic Drifts](@ref), and more.
- **Analytic Fields**: Dipole fields, magnetic bottles, tokamaks.
- **Advanced Features**: GPU acceleration, co-evolution of particles and fields, custom saving callbacks.
