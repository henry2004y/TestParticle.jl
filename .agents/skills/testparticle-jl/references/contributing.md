# Contributing workflows

`AGENTS.md` and `CONTRIBUTING.md` at the repository root are authoritative. This
file adds the details that are easy to get wrong.

## Rules that override generic habits

- Maximum line length 92 characters, enforced by Runic.jl, not JuliaFormatter.
  There is no `.JuliaFormatter.toml`.
- Explicit `;` before keyword arguments in calls.
- Avoid `Any` in function signatures.
- Comment only where the code does not already say it.
- Never edit `Manifest.toml` by hand. Use `Pkg.add` and `Pkg.rm`; `Project.toml`
  may be edited for the `[compat]` section.
- Commit titles use `component: Brief summary`, followed by a short prose body.
  Reference external work with `owner/repo#123`.

## Running tests

Full suite through the workspace test project:

```bash
julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=test test/runtests.jl
```

Equivalently, in the root environment: `julia --project=.` then `] test`.

A single test file, which is the fast loop while iterating:

```bash
julia --project=test -e 'include("test/test_canonical.jl")'
```

`test/runtests.jl` includes `test_common.jl` first, then wraps a subset of files
in `@testset "TestParticle.jl"` and includes the rest at top level.
`test/test_Makie.jl` runs only when `makie` appears in `ARGS`.

Sub-package tests are separate projects:

```bash
julia --project=lib/BorisPushers -e 'using Pkg; Pkg.test()'
julia --project=lib/TestParticleC -e 'using Pkg; Pkg.test()'
```

CI runs the main suite on Linux, macOS and Windows with
`JULIA_NUM_THREADS = 2`, plus a `test-subpackage` job for the two `lib/`
packages.

## Writing a test

Create `test/test_<feature>.jl` as a self-contained module with a `test_`
prefix, then include it from `test/runtests.jl`:

```julia
module test_feature

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using Test
using TestParticle
using ..test_common: uniform_B, max_rel_diff

@testset "feature" begin
    let
        # ...
    end
end

end # module test_feature
```

The `isdefined(Main, :test_common)` guard keeps the file runnable standalone,
and `Base.include(Main, ...)` is required because `include` inside the module
would nest `test_common` under `test_feature`. Prefer `let` blocks to avoid
leaking names, and reuse the shared field fixtures rather than defining fresh
closures, because `prepare` specialises on the field function type and new
closures trigger recompilation.

`test/test_api_surface.jl` enumerates every exported and public name. Any API
change must be reflected there, otherwise the suite fails.

## Registering documentation

- Tutorials: add a page under `docs/src/tutorial/` and register it in the
  `pages` list of `docs/make.jl:155-168`.
- Examples: add the `.jl` (or `.md`) under `docs/examples/<category>/` and add
  the filename to the matching order array in `docs/make.jl:24-90`
  (`drifts_order`, `analytic_order`, `applications_order`, `features_order`).
  Unregistered files are neither executed nor rendered.
- `docs/examples/*.jl` are Literate sources executed at build time with
  `execute = true`, and their output is cached by SHA256. Set
  `DOCS_FORCE_REBUILD=true` to regenerate everything.
- API pages are generated with `@autodocs` in `docs/src/api.md`, in two blocks:
  `Private = false` for the public API and `Public = false` for internals. There
  is no doctest setup, so `jldoctest` blocks are not verified.
- Build locally with `julia --project=docs docs/make.jl`.

## Benchmarks

`benchmark/benchmarks.jl` is the AirspeedVelocity entry point used by
`benchmark.yml`; `benchmark/run_scaling.jl` and `ScalingCommon.jl` check scaling,
and `submit_slurm.sh` runs scaling jobs on a cluster. Local ad hoc measurements
use BenchmarkTools from `benchmark/Project.toml`.

## Python bindings

`bindings/python` wraps the library through juliacall; `publish-python.yml` and
`python-wrapper.yml` drive its CI. Changes to the public Julia API may require
updating `bindings/python/src/testparticle/__init__.py`.
