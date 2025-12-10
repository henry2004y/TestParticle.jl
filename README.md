<img src="https://github.com/henry2004y/TestParticle.jl/blob/master/docs/src/assets/TP-logo.png" width="300">

---

# TestParticle.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://henry2004y.github.io/TestParticle.jl/dev)
[![Coverage](https://codecov.io/gh/henry2004y/TestParticle.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/henry2004y/TestParticle.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10149789.svg)](https://doi.org/10.5281/zenodo.10149789)

This package provides test particle tracing in an analytic or numerical electromagnetic field via [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) and native solvers that are compatible with the DifferentialEquations interface.

## Installation

In the Julia REPL,

```julia
julia> ]
pkg> add TestParticle
```

Visualization via [Makie](https://makie.juliaplots.org/stable/) and [Plots](https://docs.juliaplots.org/stable/) are supported via recipes. Please refer to each visualization library's documentation for installations.

## Usage

TestParticle.jl is designed to work together with [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).
For example, a proton in a static magnetic field can be traced in SI units via

```julia
using TestParticle, OrdinaryDiffEq, StaticArrays
# Magnetic field
B(x) = SA[0, 0, 1e-8]
# Electric field
E(x) = SA[0,0, 0.0, 0.0]
# Initial conditions
stateinit = let x0 = [1.0, 0.0, 0.0], v0 = [0.0, 1.0, 0.1]
   [x0..., v0...]
end
# Time span
tspan = (0, 20)
# Assemble particle + fields
param = prepare(E, B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
# Trace trajectory and save positions & velocities
sol = solve(prob, Vern9())
```

Native Boris particle pusher also follows a similar interface:

```julia
dt = 3e-11 # fixed time step
savestepinterval = 10
prob = TraceProblem(stateinit, tspan, param)
sol = TestParticle.solve(prob; dt, savestepinterval)[1]
```

For plotting with Makie,

```julia
using GLMakie

plot(sol, idxs=(1, 2, 3))
```

More tutorials and examples can be found in the [doc](https://henry2004y.github.io/TestParticle.jl/dev/).
