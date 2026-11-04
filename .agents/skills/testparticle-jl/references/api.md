# TestParticle.jl API reference

Signatures below are copied from the source docstrings; line references point at
the current definitions.

## Units and constants

Default is SI. Constants live in `src/utility/constants.jl` and are exported as
`public`: `qₑ`, `mₑ`, `qᵢ`, `mᵢ`, `c`, `μ₀`, `ϵ₀`, `kB`, `Rₑ`, `BMoment_Earth`,
`eV`.

Species: `Proton`, `Electron`, `Ion(m̄, q̄)`, looked up through `SpeciesDict`.
`Species{M, Q}` stores `m` and `q` (`src/types.jl:6`).

Normalized units use the same equations with `q = 1, m = 1` and `c = 1` for the
relativistic case.

## Field specification

`Field{itd, F}` (`src/prepare.jl:36`) wraps a callable and records time
dependence. Construction from a plain function dispatches on
`is_time_dependent(f) = applicable(f, zeros(3), 0.0) || applicable(f, zeros(6), 0.0)`
(`src/prepare.jl:8`).

Accepted field arguments:

- a function `f(x)` returning a 3-vector (`SVector` preferred, `Vector` works);
- a function `f(x, t)` returning a 3-vector;
- an array of size `(3, size(grid)...)`, wrapped in an interpolator;
- `nothing`, converted to `ZeroField()`;
- a `Field` or `ZeroField` passed through unchanged.

`prepare_erased` replaces field closures with `ErasedField` to cut dispatch
permutations for static compilation.

## Tracing parameters

```julia
prepare([grid...;] E = nothing, B = nothing, F = nothing, kwargs...)
    -> (; q2m, m, E, B, F[, c])                    # src/prepare.jl:155
prepare_erased([grid...;] E, B, F, kwargs...)      # src/prepare.jl:170
prepare_canonical(; phi, A, grad_phi, grad_A, kwargs...)
    -> (; q, m, c, pf)                             # src/prepare.jl:233
prepare(pf::PotentialField; species, q, m, c, type)  # src/prepare.jl:218
prepare_gc(xv, [grid...;] E, B; species, q, m, t = 0, kw...)
    -> (state_gc, (; q, q2m, μ, E, B))             # src/gc/gc.jl:49
```

`prepare` keywords: `species = Proton`, `q`, `m` (override `species`), `c`
(speed of light for relativistic equations), `type` (number type for `q2m`,
`m`, `c`, e.g. `Float32`), `order` (interpolation order 0, 1, 3), `bc`
(FastInterpolations boundary, default `FillExtrap(NaN)`; `ClampExtrap()` in 1D),
`store = StorePolicy()`, `gridtype` (`CartesianGrid`, `RectilinearGrid`,
`StructuredGrid`), `erased`.

Grid conventions: pass `x` (1D along `dir`), `x, y` (2D) or `x, y, z` (3D) as
sorted vectors, or a Meshes.jl grid. A `LogRange` radial coordinate implies
`StructuredGrid` with spherical components `(Br, Bθ, Bϕ)`. Grid vectors must be
sorted or `prepare` throws.

`PotentialField(; phi, A, grad_phi, grad_A, erased)` (`src/prepare.jl:191`)
holds the canonical potentials. Missing gradients are computed with
`ForwardDiff`; `grad_A(x, t)[i, j] = ∂A_j/∂x_i`.

## Problems and state layout

`TraceProblem{M, ...} <: AbstractODEProblem` (`src/types.jl:28`), with
`M ∈ FullOrbit | Canonical | GuidingCenter`. Aliases:
`TraceCanonicalProblem = TraceProblem{Canonical}`,
`TraceGCProblem = TraceProblem{GuidingCenter}`.

```julia
TraceProblem(u0, tspan, p; relativistic = false, prob_func = DEFAULT_PROB_FUNC)
TraceProblem(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
TraceProblem{M}(...)        # src/problem.jl:82-103
```

State layouts and parameter tuples (`src/problem.jl:3-5`):

| Model | `u0` | `p` |
| --- | --- | --- |
| `FullOrbit` | `[x..., v...]` (or `γv` when relativistic) | `(; q2m, m, E, B, F)` |
| `Canonical` | `[x..., p...]` | `(; q, m, c, pf)` |
| `GuidingCenter` | `[X..., vpar]` | `(; q, q2m, μ, E, B)` |

Positional tuples are named by `_named_params`, so documented tuple layouts
still work. Only particle constants (`q`, `q2m`, `m`, `c`, `μ`) are promoted to
the state's float type; other entries keep theirs.

Default equations (`src/problem.jl:44-63`): in-place for mutable arrays,
out-of-place for `StaticArray`s, relativistic when `relativistic = true`.
`remake` re-selects them when `u0` changes mutability.

Equations: `trace!`/`trace` (`src/equations.jl:47,54`),
`trace_relativistic!`/`trace_relativistic` (`src/equations.jl:67,74`),
`trace_canonical!`/`trace_canonical` (`src/canonical.jl:24-28`),
`trace_gc!`/`trace_gc` (`src/gc/equations.jl:100,114`), plus the drift
variants `trace_gc_drifts!`, `trace_gc_flr!`, `trace_gc_exb!`.

## Solvers

SciML algorithms (`Tsit5`, `Vern7`, `Vern9`, symplectic integrators, ...) all
work through the standard interface:

```julia
sol = solve(prob, Vern7())
sol = solve(prob, Vern7(), EnsembleThreads(); trajectories = 100, seed = 1)
```

`solve(prob::TraceProblem, alg, ensemblealg; trajectories = 1, seed = nothing,
safetycopy = false, kw...)` is defined at `src/problem.jl:117`.

Native pushers. `Boris`, `MultistepBoris{N}` and `GCRK4` need a fixed `dt`; the
adaptive variants pick their own steps, and `GCRK45` treats `dt` as an optional
initial step.

| Algorithm | Package / location | Notes |
| --- | --- | --- |
| `Boris()` | `lib/BorisPushers` | second-order, non-relativistic |
| `AdaptiveBoris(; safety)` | `lib/BorisPushers` | step from the local gyroperiod |
| `MultistepBoris4(; n)` | `lib/BorisPushers` | higher-order gyrophase correction |
| `AdaptiveMultistepBoris{4}(; n, safety)` | `lib/BorisPushers` | adaptive variant |
| `GCRK4()` | `src/gc/gc_solver.jl:8` | fixed-step RK4 for `TraceGCProblem` |
| `GCRK45(; abstol, reltol)` | `src/gc/gc_solver.jl:15` | adaptive Dormand-Prince for `TraceGCProblem` |
| `AdaptiveHybrid(; dtmax, ...)` | `src/hybrid.jl:37` | switches GC and full orbit on adiabaticity |

`Boris` types subtype `OrdinaryDiffEqAlgorithm`; `solve` for them lives in
`ext/TestParticleBorisPushersExt.jl`, so `using BorisPushers` is required.

`AdaptiveHybrid` takes the **full-orbit** state, not a guiding-center state:
`solve(TraceProblem(u0, tspan, p), AdaptiveHybrid(; dtmax))`. `dtmax` is
required; other keywords are `dtmin`, `safety_fo`, `abstol`, `reltol`,
`maxiters`, `check_interval`, `threshold` (sets both hysteresis thresholds),
`threshold_gc_to_fo`, `threshold_fo_to_gc` (must be `<=` `threshold_gc_to_fo`),
`adiabaticity ∈ (:curvature, :gradB, :both, :jacobian)`, and
`save_adiabaticity`.

Solution-side keywords handled by the native path: `dt`, `saveat`,
`save_start`, `save_end`, `save_everystep`, `save_fields`, `save_work`,
`isoutside` (with `TerminateOutside`), `maxiters`, `seed`.

## Post-processing helpers

Fields and geometry:

- `get_fields(sol)` recomputes `E` and `B` at saved states
  (`src/utility/utility.jl:434`).
- `get_magnetic_properties(x, t, Bfunc)` returns `(; B, ∇B, κ, b̂, Bmag, JB)`
  (`src/utility/field_geometry.jl:35`), with `JB[i, j] = ∂Bᵢ/∂xⱼ`.
- `get_curvature_radius(x, t, Bfunc)`, `jacobian(f, xu, t)`,
  `derivative_t(f, xu, t)`.

Energy and work:

- `get_work_rates(xv, p, t; relativistic, magnetic_properties, E)` and
  `get_work_rates(sol)` (`src/utility/work.jl:45,109`).
- `get_work_channels(xv, p, t)` returns `[P_par, P_fermi, P_grad, P_betatron]`.
- `get_energy(sol; m, q)`, `energy2velocity(Ek; m, q)`.

Gyro quantities: `get_gyrofrequency(B; q, m)`, `get_gyroperiod(B; q, m)`,
`get_gyroradius(V, B; q, m)`, `get_gyroradius(sol, t)`,
`get_adiabaticity(...)`, `adiabaticity_components(...)`.

Guiding center conversion: `get_gc`, `get_gc_velocity`, `full_to_gc`,
`gc_to_full`, `velocity_to_canonical`, `canonical_to_velocity`,
`canonical_hamiltonian`.

Boundaries and fluxes: `get_particle_crossings`, `get_first_crossing`,
`get_particle_flux`, `TerminateOutside(isoutofdomain)`.

I/O: `save_trajectory(filepath, sol; states_matrix = true)` and variants for
`(t, u)` vectors, `EnsembleSolution`, and `(; u, t)` named tuples, defined in
`ext/TestParticleJLD2Ext.jl`. Requires `using JLD2`.

Saving hooks: `orbit`, `monitor`, `save_trajectory` are declared in
`src/TestParticle.jl:89-101` for extensions and plotting recipes.

## Sampling and phase space

- `sample_unit_sphere`, `sample_maxwellian`, `sample_velocity_ball`,
  `generate_sphere`.
- Distributions: `Maxwellian`, `BiMaxwellian`, `Kappa`, `BiKappa`, forwarded to
  `VelocityDistributionFunctions.jl` through `ext/TestParticleVDFExt.jl`.
- `bin_velocity_space`, `bin_centers`, `project_vdf`, `analytic_projection`,
  `relative_l2`, `velocity_moments`, `vdf_grid_problem`, `vdf_backward`,
  `vdf_backward_trace`, `vdf_forward_trace`, `refine_vdf_window`, `embed_vdf`.

## Interpolation and grids

`build_interpolator(gridtype, f, x..., order, bc; kw...)`,
`get_dx!`, `get_dv!`, `LazyTimeInterpolator`. GPU grids:
`GPUGrid`, `GPUUniformAxis`, `GPUNonUniformAxis`. Spherical helpers:
`sph2cart`, `cart2sph`, `sph2cartvec`, `cart2sphvec`. Grid construction for
Meshes.jl grids goes through `makegrid` and `get_cell_centers` in
`ext/TestParticleMeshesExt.jl`.
