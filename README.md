# TestParticle.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://henry2004y.github.io/TestParticle.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://henry2004y.github.io/TestParticle.jl/dev)
[![Coverage](https://codecov.io/gh/henry2004y/TestParticle.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/henry2004y/TestParticle.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10149789.svg)](https://doi.org/10.5281/zenodo.10149789)

<img src="https://github.com/henry2004y/TestParticle.jl/blob/master/docs/src/assets/logo.png" align="right" style="padding-left:10px;" width="100"/>

This package provides test particle tracing in an analytical/numerical field via [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) and native solvers that are compatible with the DifferentialEquations interface.

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
sol = solve(prob, Vern7())
```

Native Boris particle pusher also follows a similar interface:

```julia
dt = 3e-11 # fixed time step
prob = TraceProblem(stateinit, tspan, param)
# Standard Boris solver
sol = TestParticle.solve(prob, Boris(); dt, savestepinterval=10)[1]
```

Besides the standard Boris method, we also support various advanced Boris solvers:

- **Adaptive Boris**: Automatic time step selection based on local gyroperiod.
- **Multistep Boris**: Sub-cycling and higher-order gyrophase correction ([Zenitani & Kato 2025](https://arxiv.org/abs/2505.02270)).

```julia
# Adaptive Boris
sol_adaptive = TestParticle.solve(prob, AdaptiveBoris(safety=0.1))[1]

# 4th-order Multistep Boris with 2 substeps
sol_multi = TestParticle.solve(prob, MultistepBoris4(n=2); dt)[1]

# Adaptive 4th-order Multistep Boris
sol_adaptive_multi = TestParticle.solve(prob, AdaptiveMultistepBoris{4}(n=2, safety=0.1))[1]
```

For plotting with Makie,

```julia
using GLMakie

plot(sol, idxs=(1, 2, 3))
```

More tutorials and examples can be found in the [documentation](https://henry2004y.github.io/TestParticle.jl/dev/).
