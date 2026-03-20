# Features Walkthrough

`TestParticle.jl` provides a range of features to simplify and enhance your particle tracing simulations. This section provides a walkthrough of these capabilities with links to specific demos.

## Unit Conversions and Scaling

By default, SI units are used. However, you can define your own unit systems by setting reference scales for mass, charge, magnetic field, length, and time.

- **Check out the demos on unit conversions**:
    - [Dimensionless Units](@ref Dimensionless-Units-and-Normalization)
    - [Unitful Integration](@ref Unit-support)

## Tracing Backwards in Time

Tracing a particle backwards is as simple as setting a time span where the end time is less than the start time, e.g., `tspan = (0.0, -1.0)`. This is particularly useful for source identification in plasma simulations.

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

## Boris Pusher

The native Boris particle pusher is optimized for performance and energy conservation in magnetic fields.

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

- **Drifts**: ExB, Gradient-B, Curvature, and more.
- **Analytic Fields**: Dipole fields, magnetic bottles, tokamaks.
- **Advanced Features**: GPU acceleration, co-evolution of particles and fields, custom saving callbacks.
