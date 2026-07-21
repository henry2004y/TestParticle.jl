# Features Walkthrough

`TestParticle.jl` provides a range of features to simplify and enhance your particle tracing simulations. This section provides a walkthrough of these capabilities with links to specific demos.

## Unit Conversions and Scaling

By default, SI units are used. However, you can define your own unit systems by setting reference scales for mass, charge, magnetic field, length, and time.

- **Check out the demos on unit conversions**:
    - [Dimensionless Units](@ref Dimensionless-Units-and-Normalization)
    - [Unitful Integration](@ref Unit-support)

## Tracing Backwards in Time

Tracing a particle backwards is as simple as setting a time span where the end time is less than the start time, e.g., `tspan = (0.0, -1.0)`. For fixed time step methods, you also need to set `dt` to be negative (e.g. standard Boris pusher). This is particularly useful for source identification in plasma simulations.

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

- [Multithreaded parallelization](@ref) (`EnsembleSerial` / `EnsembleThreads`) runs entirely inside one Julia session — no extra setup is needed, and the achievable speedup depends on `Threads.nthreads()`.
- [Distributed parallelization](@ref) (`EnsembleDistributed` / `EnsembleSplitThreads`) requires adding worker processes up front and loading the required packages on every worker, but scales beyond a single machine.

For GPU-based massively parallel tracing, see the [GPU Ensemble Tracing](@ref) example.

## [Boris Pusher](@id features-boris-pusher)

The native Boris pusher is optimized for performance and energy conservation in magnetic fields. `TestParticle.jl` provides a whole family of Boris solvers — standard, Multistep, Hyper, and Adaptive — described in [Boris Pusher](@ref Boris-Pusher).

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
