# Numerical Methods

Efficient and accurate particle tracing requires a good understanding of how to set up initial conditions and choose the right numerical algorithms.

## Setting Initial Conditions

The initial state of a particle is typically a 6-element vector or an `SVector` representing its position $\mathbf{x}$ and its velocity $\mathbf{v}$ (or a related quantity for relativistic cases).

### Particle Types

The `species` keyword of `prepare` takes a `Particle` from
[PlasmaPopulations.jl](https://github.com/JuliaPlasma/PlasmaPopulations.jl), re-exported by
`TestParticle.jl`: `Proton`, `Electron`, `Positron`, an ion or molecule from its formula and
charge number, or any charge and mass.

```julia
Particle(:O, 1)                          # O⁺, standard atomic weight
Particle("{4}He", 2)                     # α particle
Particle(; q = -1.0e-16, m = 1.0e-15)    # charged grain [C], [kg]
```

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
    With `TraceProblem{RelativisticOrbit}(stateinit, tspan, param)` (or `trace_relativistic!`), the solver expects the last three elements of the state to be $\gamma \mathbf{v}$.

## Choice of numerical algorithms

Explicit Runge--Kutta solvers (such as `Tsit5` or `Vern7`/`Vern9`) do not preserve phase-space volume and accumulate secular energy drift over hundreds of gyroperiods unless tight tolerances or small steps are enforced. Prefer volume-preserving pushers from [BorisPushers.jl](https://github.com/henry2004y/TestParticle.jl/tree/master/lib/BorisPushers) for long full-orbit integrations, and reserve SciML ODE solvers for guiding-center, canonical, custom-callback, or external-force problems:

| Plasma Regime | Recommended Solver(s) | Notes |
| :--- | :--- | :--- |
| **Non-relativistic full orbit** | `MultistepBoris6(; n = 4)` (phase accuracy), `Boris()` / `AdaptiveBoris()` (throughput / strong $\nabla B$) | Conserves energy to round-off; leapfrog pushers carry a half-step velocity offset in explicitly time-varying $\mathbf{E}(t)$. |
| **Relativistic full orbit ($\mathbf{E} \neq 0$)** | `HigueraCary()` or `AdaptiveHigueraCary()` (or `Vay()` / `AdaptiveVay()`) | Exact $\mathbf{E}\times\mathbf{B}$ force cancellation; `HigueraCary` is also volume-preserving. |
| **Relativistic full orbit ($\mathbf{E} = 0$)** | `MultistepBoris6(; n = 4)` or `Boris()` | High gyrophase fidelity in pure-$\mathbf{B}$ turbulence. |
| **Globally adiabatic ($\rho_L \ll L_B$)** | `TraceGCProblem` with `Vern7()` | Averages out gyromotion; currently non-relativistic and omits the $\partial_t\hat{\mathbf{b}}$ geometric drift. |
| **Multi-scale (adiabatic + non-adiabatic)** | `AdaptiveHybrid(; threshold = 0.1, dtmax)` on `TraceProblem` | Automatically switches between `AdaptiveBoris` (FO) and `Vern7()` (GC). |
| **Canonical / generic force $\mathbf{F}$ / continuous callbacks** | `Vern7()`, `Vern9()`, or symplectic/implicit SciML ODE (`ImplicitMidpoint()`) | Supports `TraceCanonicalProblem`, non-EM forces $\mathbf{F}$, and arbitrary SciML callbacks. |
| **Large GPU ensembles ($N \gtrsim 10^5$)** | `Boris()`, `MultistepBoris*`, `Vay()`, `HigueraCary()` with `raw_output = true` | Pass a contiguous $N \times 6$ `u0` matrix to avoid per-particle host wrapper overhead. |

By default `DifferentialEquations.jl` applies `Tsit5` to an ODE problem. If only `OrdinaryDiffEqVerner.jl` (or `OrdinaryDiffEq.jl`) is imported, pass the solver explicitly (e.g., `Vern7()` or `Vern9()`, tuning `reltol`/`abstol` as needed).

## Native Boris solver

You can also try out the native implementation of the Boris method, provided by [BorisPushers.jl](https://github.com/henry2004y/TestParticle.jl/tree/master/lib/BorisPushers). The Boris pusher is specifically designed for particle tracing in magnetic fields and has the property of exactly conserving energy in a static, uniform magnetic field. Once `BorisPushers` is loaded, it follows a similar interface for ease of adoption:

```julia
using BorisPushers
dt = 1e-2 # fixed time step, about 1/650 of the gyroperiod (≈ 6.6 s)
prob = TraceProblem(stateinit, tspan, param)
sol = solve(prob, Boris(); dt, saveat = 10 * dt)
```

For the full family of Boris solvers (Multistep, Hyper, and Adaptive), see [Boris Pusher](@ref Boris-Pusher). When you want the accuracy of that family without tuning it yourself, `MultistepBoris6(; n = 4)` is the recommended default.

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
