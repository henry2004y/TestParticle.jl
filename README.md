# TestParticle.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://henry2004y.github.io/TestParticle.jl/dev)
[![Coverage](https://codecov.io/gh/henry2004y/TestParticle.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/henry2004y/TestParticle.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10149789.svg)](https://doi.org/10.5281/zenodo.10149789)

This package provides test particle tracing in an analytic or numerical electromagnetic field via [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).

## Installation

In the Julia REPL,

```julia
julia> ]
pkg> add TestParticle
```

Visualization via [Makie](https://makie.juliaplots.org/stable/), [Plots](https://docs.juliaplots.org/stable/), and [PyPlot](https://github.com/JuliaPy/PyPlot.jl) are supported. Please refer to each visualization library's documentation for installations.

## Usage

TestParticle.jl is designed to work together with [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl). Native particle pusher also follows a similar interface.

For example, a proton in a static magnetic field can be traced via

```julia
using TestParticle, OrdinaryDiffEq, StaticArrays
# Magnetic field
B(x) = SA[0, 0, 1e-8]
# Electric field
E(x) = SA[0,0, 0.0, 0.0]
# Initial conditions
x0 = [1.0, 0.0, 0.0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 20)
# Assemble particle + fields
param = prepare(E, B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
# Trace trajectory and save positions & velocities
sol = solve(prob, Vern9())
```

For plotting with Makie,

```julia
using GLMakie

plot(sol)
```

More tutorials and examples can be found in the [doc](https://henry2004y.github.io/TestParticle.jl/dev/).
