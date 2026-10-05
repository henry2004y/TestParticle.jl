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
| `MultistepBoris2(; n=1)` | Multi-step Boris with `n` sub-cycles (Zenitani & Kato 2025) |
| `MultistepBoris4(; n=1)` | 4th-order Hyper-Boris |
| `MultistepBoris6(; n=1)` | 6th-order Hyper-Boris |
| `AdaptiveBoris(; safety=0.1)` | Boris with gyroperiod-based adaptive time step |
| `AdaptiveMultistepBoris{N}(; n=1, safety=0.1)` | Multi-step Boris with adaptive time step |

All solvers subtype `AbstractBoris` and are `isbits` structs, making them safe
to pass into GPU kernels. A stateless step API (`advance_boris`, `update_velocity_half`,
`update_velocity_node`) is provided for use inside
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) kernels.

## Relation to TestParticle.jl

`BorisPushers` owns the physics and the SciML algorithm interface.
[TestParticle.jl](https://github.com/henry2004y/TestParticle.jl) depends on it
and wraps it with particle-tracing-specific features: the `TraceProblem`
container, boundary callbacks, and field/work output columns.

The package can be used independently of `TestParticle.jl` for any problem that
can be expressed as an `ODEProblem` or `EnsembleProblem`. The parameter container
only needs to satisfy three accessor functions: `get_q2m`, `get_EField`, and `get_BField`.

```julia
using BorisPushers, StaticArrays

q2m, m = 1.0, 1.0
Efunc(x, t) = SA[0.0, 0.0, 0.0]
Bfunc(x, t) = SA[0.0, 0.0, 1.0]

param = (q2m, m, Efunc, Bfunc)
u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
tspan = (0.0, 10.0)

prob = ODEProblem((u, p, t) -> nothing, u0, tspan, param)
sol = solve(prob, Boris(); dt = 0.1)
```

## GPU & Hardware Accelerator Support

`BorisPushers` provides `EnsembleKernel` for device execution across CPU and GPU
architectures (CUDA, AMDGPU, oneAPI, Metal, or multi-threaded CPU) via
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl):

```julia
using BorisPushers, StaticArrays, KernelAbstractions

# Run on CPU threads or pass a GPU backend like CUDA.CUDABackend()
ens = EnsembleKernel(CPU())
sol = solve(prob, Boris(), ens; dt = 0.1, trajectories = 1000)
```

## Documentation

Full documentation, including algorithm descriptions, solver options,
ensemble usage, and GPU ensemble tracing, is part of the
[TestParticle.jl documentation](https://henry2004y.github.io/TestParticle.jl/dev/tutorial/boris.html).

## License

BorisPushers is released under the MIT License. See [LICENSE](LICENSE).
