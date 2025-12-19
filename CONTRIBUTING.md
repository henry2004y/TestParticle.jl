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

We use [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) with the **SciML** style and **3-space indentation**.

To format your code:

```julia
using JuliaFormatter
format(".")
```

## Building Documentation

To build the documentation locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated HTML files will be in `docs/build/`.

---

Happy coding!
