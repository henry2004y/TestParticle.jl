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

There are two primary ways to trace multiple particles simultaneously:

1. Extracting the solution in a loop with varying initial conditions. See the example [Ensemble Tracing](@ref).
2. Constructing the [Ensemble Simulations](https://diffeq.sciml.ai/stable/features/ensemble/). One example can be found [here](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_ensemble.jl). However, note that by default the ensemble type replicates the parameters for each solution, which is very memory inefficient for tracing in a numeric field. 
We need to set `safetycopy=false` to make the field as a reference in the parameter of each trajectory.

- **Example**: [Ensemble Tracing](@ref Ensemble-Tracing)

!!! warning
    When defining a `prob_func` for an `EnsembleProblem`, avoid in-place mutation of the `ODEProblem` (e.g., `prob.u0 = ...`), especially when using `StaticArrays`. SciMLBase v2.0+ may issue a warning about "Mutation of ODEProblem detected". Instead, use `remake(prob; u0=...)` to create a new problem instance with updated parameters.

## Boris Method

The native Boris particle pusher is optimized for performance and energy conservation in magnetic fields.

- **Example**: [Boris Method](@ref Boris-Method)

## Field Interpolation

`TestParticle.jl` supports tracing in numerical fields defined on various grids (Cartesian, Rectilinear, Structured).

- **Example**: [Field Interpolation](@ref Field-Interpolation)

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
