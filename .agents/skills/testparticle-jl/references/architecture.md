# Architecture

## Load order

`src/TestParticle.jl` is the only place that defines what the package contains.
It imports dependencies, lists `export` and `public` names, and then includes
files in a fixed order that the rest of the source relies on:

```
types.jl
utility/utility.jl -> utility/{constants,zero,phasespace}.jl
utility/field_geometry.jl
utility/work.jl
utility/interpolation.jl
utility/gpu_grid.jl
sampler.jl
prepare.jl -> erased.jl
saveat.jl
gc/{gc,gc_solver,rk4_gc,rk45_gc,equations}.jl
equations.jl
canonical.jl
problem.jl
hybrid.jl
deprecated.jl
precompile.jl
```

Adding a new file means adding an `include` here, at a position after the
definitions it uses.

## Type map

| Type | File | Role |
| --- | --- | --- |
| `AbstractField{itd}` | `src/types.jl:1` | field callable supertype, `itd` is time-dependence |
| `Field{itd, F}` | `src/prepare.jl:36` | wraps a function or interpolator |
| `ErasedField`, `ErasedFunction` | `src/erased.jl` | type-erased wrappers for static compile |
| `ZeroField`, `ZeroVector` | `src/utility/zero.jl` | zero-field fast paths |
| `Species{M, Q}` | `src/types.jl:6` | mass and charge |
| `FullOrbit`, `Canonical`, `GuidingCenter` | `src/types.jl:24-26` | model tags |
| `TraceProblem{M, ...}` | `src/types.jl:28` | `AbstractODEProblem` plus `prob_func` |
| `PotentialField` | `src/types.jl:68` | canonical potentials |
| `CartesianGrid`, `RectilinearGrid`, `StructuredGrid` | `src/types.jl:76-78` | abstract stubs bound by the Meshes extension |

`TraceProblem` fakes the `kwargs` property that SciMLBase solvers read, through
`Base.propertynames` and `Base.getproperty` (`src/types.jl:55-60`).

## Subprojects in `lib/`

Neither subproject is a workspace member; the root `Project.toml` lists only
`["test", "docs"]`. Both are wired in through path sources in
`test/Project.toml`.

- `lib/BorisPushers` (package `BorisPushers`, v0.2.0): native Boris family as
  `OrdinaryDiffEqAlgorithm` subtypes.
  - `algorithms.jl`: `AbstractBoris`, `Boris`, `MultistepBoris{N}` (2, 4, 6),
    `AdaptiveBoris`, `AdaptiveMultistepBoris{N}`.
  - `boris_caches.jl`, `boris_perform_step.jl`, `boris_controller.jl`: the
    OrdinaryDiffEqCore integration points.
  - `boris_solve.jl`: fast `solve!` path.
  - `kernels.jl`, `ensemble_kernel.jl`: KernelAbstractions GPU ensemble.
- `lib/TestParticleC` (package `TestParticleC`): dependency-light tracing kernel
  callable from C, built from StaticArrays and BorisPushers only.

## Extension mechanism

Stub functions are declared in the main module and extended in `ext/`:

- `function orbit end`, `function monitor end`, `function save_trajectory end`
  (`src/TestParticle.jl:89-93`). `save_trajectory` has a fallback method that
  throws an error telling the user to load JLD2.
- The `solve` methods for Boris algorithms and `AdaptiveHybrid` live in
  `ext/TestParticleBorisPushersExt.jl`. That file converts a
  `TraceProblem{FullOrbit}` into a fully specialised `ODEProblem`, gives it a
  no-op right-hand side, and rebuilds fields and work rates afterwards.

Weak dependencies and extension names are declared in the root `Project.toml`
under `[weakdeps]` and `[extensions]`. There are five:

| Extension | Weak dep | Provides |
| --- | --- | --- |
| `TestParticleBorisPushersExt` | `BorisPushers` | `solve` for Boris algorithms and `AdaptiveHybrid` |
| `TestParticleJLD2Ext` | `JLD2` | `save_trajectory` methods |
| `TestParticleMeshesExt` | `Meshes` | grid types, `makegrid`, `build_interpolator` dispatch |
| `TestParticleRecursiveArrayToolsExt` | `RecursiveArrayTools` | `trace`/`trace!` for `ArrayPartition` |
| `TestParticleVDFExt` | `VelocityDistributionFunctions` | `sample_maxwellian`, `Maxwellian`, `Kappa`, ... |

## Recipes

### Add a post-processing helper

1. Put the implementation next to related code (`utility/utility.jl` for
   diagnostics, `utility/work.jl` for energetics, `utility/field_geometry.jl`
   for curvature and gradient quantities).
2. Add the name to the `public` list in `src/TestParticle.jl:44-67`, or to
   `export` only if it is a first-class user API.
3. If it has methods taking a solution, also accept `AbstractODESolution` and
   read `sol.prob.p`, following `get_work_rates(sol)` and
   `get_adiabaticity(sol)`.
4. Add the name to `test/test_api_surface.jl`; that test enumerates every
   exported and public name.

### Add a solver

Native solvers are dispatched from `solve(prob::TraceProblem, alg, ...)`.
For a SciML-style algorithm, subtype `OrdinaryDiffEqAlgorithm` in a subproject
and bridge it in an extension, mirroring `BorisPushers`. For a self-contained
integrator, follow the GC path: `src/gc/gc_solver.jl` dispatches on the
algorithm, `rk4_gc.jl` and `rk45_gc.jl` hold the update kernels, and
`saveat.jl` handles output.

### Add a field type or field behaviour

Fields flow through `prepare_field` (`src/prepare.jl:91-102`). Add a method
there rather than branching inside `prepare`. Keep the `AbstractField{itd}`
parameter accurate; `is_time_dependent(::AbstractField{itd})` returns `itd`.

### Add an extension

Add the weak dependency to `[weakdeps]`, the extension entry to
`[extensions]`, and the file `ext/TestParticle<Name>Ext.jl` with a
`module TestParticle<Name>Ext` that imports the parent package explicitly.

## GPU and AD considerations

- GPU paths use KernelAbstractions and live in `lib/BorisPushers/src/kernels.jl`,
  `ensemble_kernel.jl`, and `src/utility/gpu_grid.jl`.
- Differentiability uses `ForwardDiff` for `jacobian` and `derivative_t`.
- `erased` fields and `@setup_workload`/`@compile_workload` in
  `src/precompile.jl` exist to control compile time; keep them in sync with new
  entry points.
