---
name: testparticle-jl
description: This skill should be used for any task involving TestParticle.jl, the Julia test-particle tracing package (full-orbit, relativistic, canonical, and guiding-center tracing in analytic or gridded fields). Triggers include writing or debugging particle simulations with prepare/TraceProblem/solve, choosing between SciML ODE solvers, Boris pushers, GCRK4/GCRK45 and AdaptiveHybrid, handling SI versus normalized units, and contributing code, tests, benchmarks, or docs to the TestParticle.jl repository.
---

# TestParticle.jl

## Purpose

Set up, run, and extend test-particle tracing with TestParticle.jl. The package
integrates charged-particle motion in electromagnetic fields through the SciML
`ODEProblem`/`solve` interface, adds native Boris and guiding-center solvers with
the same interface, and provides sampling, post-processing, and visualization
helpers.

Use this skill both when *using* the package in a simulation script and when
*modifying* the package in this repository.

## Workflow for a tracing task

Follow these five steps in order. Do not skip step 2, since the model choice
determines the state layout and the available solvers.

1. **Describe the fields.** A field is either a function `f(x)` or `f(x, t)`
   returning a 3-vector (`SVector` preferred), or an array of size
   `(3, size(grid)...)` interpolated on a sorted grid. Omitted fields are zero,
   so `prepare(; B)` is valid. See `references/api.md` for grid and
   interpolation keywords.
2. **Choose the model.**
   | Model | Setup | Native solvers |
   | --- | --- | --- |
   | `FullOrbit` (default) | `prepare(; E, B[, F])` | `Boris`, `MultistepBoris4`, `AdaptiveBoris`, `AdaptiveMultistepBoris{4}` |
   | Relativistic full orbit | `prepare(; E, B, c)` plus `relativistic = true` | SciML adaptive solvers only (`Boris` is non-relativistic) |
   | `Canonical` | `prepare_canonical(; phi, A)` | SciML adaptive solvers only |
   | `GuidingCenter` | `prepare_gc(xv; E, B)` | `GCRK4`, `GCRK45` |
   | GC/full-orbit hybrid | full-orbit `p` and `u0` | `AdaptiveHybrid(; dtmax)` |
3. **Build the problem.** `TraceProblem(u0, tspan, p)` with `u0 = [x..., v...]`
   for full orbit, `[x..., p...]` for canonical, `[X..., vpar]` for guiding
   center. Pass `relativistic = true` to get the relativistic equations.
4. **Solve.** `solve(prob, alg; dt, saveat, ...)`. SciML algorithms and the
   adaptive native variants choose their own steps; `Boris`, `MultistepBoris*`
   and `GCRK4` require a `dt` that resolves the local gyroperiod. Ensembles go
   through `solve(prob, alg, EnsembleThreads(); trajectories)`.
5. **Post-process.** `get_fields`, `get_work_rates`, `get_gyroradius`,
   `get_gyroperiod`, `get_particle_flux`, `get_first_crossing`, `orbit`,
   `monitor`, `save_trajectory` (needs `using JLD2`), and plotting via the
   Makie/Plots recipes.

Minimal end-to-end example (SI units, proton in a 10 nT field):

```julia
using TestParticle, OrdinaryDiffEq, StaticArrays
B(x) = SA[0.0, 0.0, 1.0e-8]
param = prepare(; B, species = Proton)
prob = TraceProblem([1.0, 0.0, 0.0, 0.0, 1.0, 0.1], (0.0, 20.0), param)
sol = solve(prob, Vern7())
```

With a native pusher (`dt` must resolve the gyroperiod, roughly 6.6 s here):

```julia
using BorisPushers
sol = solve(prob, Boris(); dt = 1e-2, saveat = 10 * dt)
```

## Units

- Default is SI throughout: `q`, `m`, `c` come from CODATA constants and field
  arrays are in tesla and volt per metre.
- Normalized units need no separate equations. Pass `q = 1, m = 1`, and
  `c = 1` for the relativistic case; see `docs/examples/features/demo_dimensionless.jl`.
- `trace_relativistic!` stores `γv`, not `v`. Convert back when comparing with
  non-relativistic output.

## Common pitfalls

- Do not pass zero fields explicitly unless the signature requires it;
  `ZeroField` is the cleaner default and takes a fast path.
- Field functions are inspected with `applicable(f, zeros(3), 0.0)` to decide
  whether they are time dependent. A function whose second argument has a
  default value silently becomes time dependent.
- Mutable `Vector` states select the in-place equations, `StaticArray` states
  the out-of-place ones. `remake` re-picks the default equations accordingly.
- Grid vectors must be sorted; array fields must have the component axis first.
- `AdaptiveHybrid` requires `using BorisPushers`, since its `solve` method lives
  in an extension.
- `save_trajectory` throws until `using JLD2` loads its extension.
- Names in the `public` list (constants such as `qᵢ` and `mᵢ`, diagnostic
  helpers such as `get_gyroperiod`, distributions, interpolation policies) are
  not exported. Import them explicitly, e.g.
  `using TestParticle: get_gyroperiod, qᵢ`.

## Repository map

- `src/` main package; `src/TestParticle.jl` is the single source of truth for
  the include order and the `export`/`public` lists.
- `src/gc/` guiding center, `src/utility/` constants, interpolation, geometry,
  work rates, phase space.
- `lib/BorisPushers/` native Boris solvers as an `OrdinaryDiffEqAlgorithm`
  subtype, `lib/TestParticleC/` dependency-light C-callable kernel. Both are
  separate packages wired into the test environment by path, not workspace
  members.
- `ext/` five weak-dependency extensions (`BorisPushers`, `JLD2`, `Meshes`,
  `RecursiveArrayTools`, `VelocityDistributionFunctions`).
- `test/` independent `module test_*` files, `docs/` Documenter + VitePress with
  Literate examples, `benchmark/` AirspeedVelocity plus BenchmarkTools,
  `bindings/python/` juliacall wrapper.

## Contributing

Read `AGENTS.md` at the repository root first; it holds the formatting and test
rules that override generic habits. Then load:

- `references/architecture.md` for where types, equations, and extensions live,
  and for the recipes "add a field type", "add a solver", "add an extension".
- `references/api.md` for signatures, keyword arguments, and helper functions.
- `references/contributing.md` for running a single test, the API-surface test,
  Runic formatting, CI jobs, docs page registration, and benchmarks.

Run `julia --project=test .agents/skills/testparticle-jl/scripts/smoke_check.jl`
from the repository root to confirm the environment still traces a proton
correctly before and after a change.
