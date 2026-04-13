# Numerical Methods

Efficient and accurate particle tracing requires a good understanding of how to set up initial conditions and choose the right numerical algorithms.

## Setting Initial Conditions

The initial state of a particle is typically a 6-element vector or an `SVector` representing its position $\mathbf{x}$ and its velocity $\mathbf{v}$ (or a related quantity for relativistic cases).

### Particle Types

`TestParticle.jl` provides several predefined species that can be passed to the `prepare` function using the `species` keyword argument:

- `Proton`: $m \approx 1.67 \times 10^{-27}$ kg, $q \approx 1.60 \times 10^{-19}$ C
- `Electron`: $m \approx 9.11 \times 10^{-31}$ kg, $q \approx -1.60 \times 10^{-19}$ C
- `Ion`: A generic ion type.

You can also define custom species using the `Species(m, q)` struct.

### Non-relativistic Case

For non-relativistic simulations, the initial state is `[x, y, z, vx, vy, vz]`.

```julia
x0 = [1.0, 0.0, 0.0]
v0 = [100.0, 0.0, 0.0]
stateinit = [x0..., v0...]
```

### Relativistic Case

In relativistic simulations, the state is typically represented by its position $\mathbf{x}$ and the quantity $\gamma \mathbf{v}$, where $\gamma$ is the Lorentz factor. `TestParticle.jl` handles the conversion to velocity internally using the `get_relativistic_v` function.

```julia
x0 = [1.0, 0.0, 0.0]
v0 = [0.8c, 0.0, 0.0] # 80% of speed of light
# For relativistic equations, the state is (x, \gamma v)
gamma = 1 / sqrt(1 - (norm(v0)/c)^2)
stateinit = [x0..., (gamma .* v0)...]
```

!!! note "Relativistic velocity"
    When using `trace_relativistic!`, the solver expects the last three elements of the state to be $\gamma \mathbf{v}$.

## Choice of numerical algorithms

By default `DifferentialEquations.jl` applies `Tsit5` to an ODE problem. If only `OrdinaryDiffEq.jl` is imported, then we need to explicitly select a solver.

However, not all solvers are guaranteed to work. For example, some cases with strong magnetic fields are tested to work only with fixed timestep algorithms like `Euler` and the Adams-Bashforth family.

Currently we recommend `Vern9` as a starting point for adaptive timestepping, with additional fine tuning by changing `reltol` if needed.

## Native Boris solver

You can also try out the native implementation of the Boris method in `TestParticle.jl`. The Boris pusher is specifically designed for particle tracing in magnetic fields and has the property of exactly conserving energy in a static, uniform magnetic field. It follows a similar interface for ease of adoption:

```julia
dt = 3e-11 # fixed time step
savestepinterval = 10
prob = TraceProblem(stateinit, tspan, param)
sol = TestParticle.solve(prob; dt, savestepinterval)[1]
```

## Boundary check

Boundary check is performed via the `callback` keyword argument or `isoutside` keyword argument for the Boris solvers. By default, no boundary check is performed, but the tracing may stop if NaN values are encountered when interpolating the EM fields or numerical instability occurs.

!!! tip "Boundary Check"
    It is recommended to use the `TerminateOutside` helper for common boundary conditions. For more complex logic, you can use `DiscreteCallback` directly.

```julia
isoutside(u, p, t) = norm(u) < Rₑ
callback = TerminateOutside(isoutside)
sol = solve(prob; callback)
```

The `callback` system is supported by both native and SciML-based solvers in `TestParticle.jl`.

!!! warning "Deprecated isoutofdomain"
    The `isoutofdomain` usage in the Boris solvers is deprecated and not recommended for the SciML-based solvers. Please use `TerminateOutside` instead.

## Solution Interpolations

All the tracing solutions come with an interpolation method to make them continuous. Depending on the order of the scheme, different orders of interpolations are applied.

Note that if the scheme comes with a "lazy" interpolation (e.g. Vern family), you can either output the full solution (default) and interpolate later, or turn off the lazy interpolation and select part of the solution to save.

```julia
# Default: lazy interpolation enabled
sol = solve(prob, Vern9())

# Selective saving with lazy interpolation disabled
sol = solve(prob, Vern9(lazy=false), save_idxs=[1,2,3])
```

The interpolated solution would be wrong if `lazy=true` and `save_idxs` is not the full state vector `[1,2,3,4,5,6]`.
