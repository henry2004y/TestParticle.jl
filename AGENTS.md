See [CONTRIBUTING.md](./CONTRIBUTING.md) for Runic formatting (92-char limit), API visibility (`export` vs `public`), and environment setup.

# Code and commits
- Avoid `Any` in function signatures, and use an explicit `;` before keyword arguments in calls.
- Add comments or docstrings only when behavior is not obvious from names and types.
- Commit titles follow `component: Brief summary` with a brief prose body; link external issues/PRs as `owner/repo#123`.
- Respect manual user edits to code you wrote; explain concerns and ask before reintroducing modified or deleted code.

# Dependencies
- Manage dependencies with `Pkg.add` and `Pkg.rm`; edit `Project.toml` directly only for `[compat]`.
- Subpackages `lib/BorisPushers` and `lib/TestParticleC` are wired into `test` and `docs` via `[sources]` in their `Project.toml` files.

# Testing and doctests
- Single test file (fast loop): `julia --project=test -e 'include("test/test_canonical.jl")'`
- Full test suite: `julia --project=test test/runtests.jl`
- Subpackages: `julia --project=lib/BorisPushers -e 'using Pkg; Pkg.test()'` and `julia --project=lib/TestParticleC -e 'using Pkg; Pkg.test()'`
- Doctests: `julia --project=docs docs/make.jl doctest`
- Write new tests as standalone `module test_<name>` files included from [`test/runtests.jl`](./test/runtests.jl), placing imports and `struct` definitions at module top level and isolating cases in `let` blocks.
- Update [`test/test_api_surface.jl`](./test/test_api_surface.jl) whenever `export` or `public` names change.
