# Contributing to TestParticle.jl

Thank you for your interest in contributing! We welcome all contributions, from bug fixes and documentation improvements to new features.

## How to Contribute

1. **Report Issues**: If you find a bug or have a feature request, please open an issue on GitHub.
2. **Submit Pull Requests**:
   - Fork the repository and create a new branch.
   - Make your changes.
   - Ensure tests pass.
   - Format your code.
   - Open a Pull Request (PR) with a clear description.

## Development Setup

To set up the development environment:

```julia
julia> ]
pkg> dev .
```

## Running Tests

To run the test suite, use the `test` project environment:

```bash
julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=test test/runtests.jl
```

## Code Style

We use [Runic.jl](https://github.com/fredrikekre/Runic.jl), with the recommended line character limit of 92.

## API Visibility

Choose visibility by role, consistently across related API families:

- **Export** tracing setup and execution: preparation functions, problem constructors,
  tracing right-hand sides (including canonical and guiding-center variants), solver
  algorithms, particle species constructors, field constructors, and `build_interpolator`.
  Also export end-to-end VDF tracing, the main `get_fields`/`get_work` result
  accessors, `orbit`/`monitor`, and the SciML ensemble/callback interfaces used
  in the tracing workflow.
- **Public but not exported** diagnostics, state and coordinate conversions, sampling
  functions and distribution types, physical constants, interpolation policies and storage
  types, and lower-level field/reconstruction building blocks. Use qualified names or
  explicit imports for these APIs.
- **Internal** implementation details that are not intended as supported user interfaces.
  Do not declare these `public` merely because tests or other package code use them.

Both exported and public-only names are supported API and carry compatibility commitments.
Moving an exported name to public-only can break code that relied on `using TestParticle`;
update examples and callers with explicit imports, and qualify documentation cross-references.

## Building Documentation

To build the documentation locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated HTML files will be in `docs/build/`.

---

Happy coding!
