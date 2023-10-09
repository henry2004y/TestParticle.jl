# ---
# title: Ensemble tracing
# id: demo_ensemble
# date: 2023-04-20
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Tracing multiple charged particles in a static EM field
# ---

# This example demonstrates tracing multiple electrons in an analytic EM field and how to
# take advantage of the multithreading support in the ODE solver. A multiproc version is
# also available. Check the official documentation of DifferentialEquations.jl for details.
# Note: as of OrdinaryDiffEq v6.31.2, parameters of the solutions in the EnsembleProblem are
# replicated for each particle solution, which is highly memory inefficient especially when
# numerical EM fields are used!

import DisplayAs # hide

using TestParticle
using OrdinaryDiffEq
using StaticArrays
using TestParticleMakie
using CairoMakie
CairoMakie.activate!(type = "png")

"Initial state perturbation for EnsembleProblem."
function prob_func(prob, i, repeat)
   remake(prob, u0=rand()*prob.u0)
end

## Initialization

function B(x)
   return SA[0, 0, 1e-11]
end

function E(x)
   return SA[0, 0, 5e-13]
end

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param = prepare(E, B, species=Electron)
tspan = (0.0, 10.0)

trajectories = 10

## Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sols = solve(ensemble_prob, Tsit5(), EnsembleThreads();
   trajectories=trajectories, save_idxs=[1,2,3])

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Electron trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data,
)

for i in eachindex(sols)
   lines!(ax, sols[i], label="$i")
end

f = DisplayAs.PNG(f) # hide