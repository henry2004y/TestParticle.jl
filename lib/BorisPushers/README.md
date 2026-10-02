# BorisPushers.jl

**BorisPushers** is a Julia sub-package of
[TestParticle.jl](https://github.com/henry2004y/TestParticle.jl) that provides
a family of Boris particle-pusher algorithms compatible with the
[SciML](https://sciml.ai/) `solve` interface.

## Purpose

The package implements several Boris-type integrators for tracing charged
particles in electromagnetic fields:

| Solver | Description |
|---|---|
| `Boris()` | Standard second-order Boris pusher |
| `MultistepBoris2(n)` | Multi-step Boris with `n` sub-cycles |
| `MultistepBoris4(n)` | 4th-order Hyper-Boris (Zenitani & Kato 2025) |
| `MultistepBoris6(n)` | 6th-order Hyper-Boris |
| `AdaptiveBoris(safety)` | Boris with gyroperiod-based adaptive time step |

All solvers are `isbits` structs, making them safe to pass into GPU kernels.
A stateless step API (`advance_boris`, `update_velocity_half`,
`update_velocity_node`) is provided for use inside
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) kernels.

## Relation to TestParticle.jl

`BorisPushers` owns the physics and the SciML algorithm interface.
[TestParticle.jl](https://github.com/henry2004y/TestParticle.jl) depends on it
and wraps it with particle-tracing-specific features: the `TraceProblem`
container, boundary callbacks, and field/work output columns.

The package can be used independently of `TestParticle.jl` for any problem that
can be expressed as an `ODEProblem`. The parameter object only needs to satisfy
three accessor functions: `get_q2m`, `get_EField`, and `get_BField`.

```julia
using BorisPushers, StaticArrays

# Build a parameter container that provides the three accessors,
# or use the default (q2m, m, Efunc, Bfunc) tuple layout.
param = (q2m, m, Efunc, Bfunc)
prob = ODEProblem((u, p, t) -> nothing, SA[x0..., v0...], tspan, param)
sol = solve(prob, Boris(); dt)
```

## Documentation

Full documentation, including algorithm descriptions, solver options,
ensemble usage, and GPU ensemble tracing, is part of the
[TestParticle.jl documentation](https://henry2004y.github.io/TestParticle.jl/dev/tutorial/boris.html).

## License

BorisPushers is released under the MIT License. See [LICENSE](LICENSE).
