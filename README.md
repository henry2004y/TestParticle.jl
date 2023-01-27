# TestParticle.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://henry2004y.github.io/TestParticle.jl/dev)
[![Build Status](https://github.com/henry2004y/TestParticle.jl/workflows/CI/badge.svg)](https://github.com/henry2004y/TestParticle.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![Coverage](https://codecov.io/gh/henry2004y/TestParticle.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/henry2004y/TestParticle.jl)

This package provides test particle tracing in an analytic or numerical electromagnetic field via [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).

## Installation

In the Julia REPL,

```julia
julia> ]
pkg> add TestParticle
```

Visualization via [Makie](https://makie.juliaplots.org/stable/), [Plots](https://docs.juliaplots.org/stable/), and [PyPlot](https://github.com/JuliaPy/PyPlot.jl) are supported. Please refer to each visualization library's documentation for installations.

## Usage

TestParticle.jl is designed to work together with [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

A proton trajectory in a static magnetic field can be traced via

```julia
using TestParticle, OrdinaryDiffEq, StaticArrays
# Magnetic field
B(x) = SA[0, 0, 1e-8]
# Electric field
E(x) = SA[0, 0, 0]
# Initial conditions
x0 = [1.0, 0.0, 0.0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 20)
# Assemble particle + fields
param = prepare(E, B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
# Trace trajectory and save positions & velocities
sol = solve(prob, Tsit5(); save_idxs=[1,2,3,4,5,6])
```

For plotting with Makie,

```julia
using TestParticleMakie, GLMakie

plot(sol)
```

More tutorials and examples can be found in [the doc](https://henry2004y.github.io/TestParticle.jl/dev/).
