# ---
# title: Ensemble tracing
# id: demo_ensemble
# date: 2024-01-23
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.0
# description: Tracing multiple charged particles in a static EM field
# ---

# This example demonstrates tracing multiple electrons in an analytic EM field and how to
# take advantage of the multithreading support in the ODE solver. A multiproc version is
# also available. Check the official documentation of DifferentialEquations.jl for details.
# In performing test particle tracing, we want to share the field information for all
# particles. This can be achieved in the ensemble problem with `safetycopy=false`.

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

"Set initial state for EnsembleProblem."
function prob_func(prob, i, repeat)
   ## If u0 is immutable (e.g. StaticArrays), use remake; otherwise we can directly modify u0
   ##remake(prob, u0=[i/3, 0.0, 0.0])
   prob.u0[4] = i / 3 

   prob
end

## Initialization

B(x) = SA[0, 0, 1e-11]
E(x) = SA[0, 0, 1e-13]

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param = prepare(E, B, species=Electron)
tspan = (0.0, 10.0)

trajectories = 3

## Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
sols = solve(ensemble_prob, Tsit5(), EnsembleThreads(); trajectories)

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
   lines!(ax, sols[i], idxs=(1,2,3), label="$i")
   ##TODO: wait for https://github.com/MakieOrg/Makie.jl/issues/3623 to be fixed!
   ax.scene.plots[9+2*i-1].color = Makie.wong_colors()[i]
end

f = DisplayAs.PNG(f) #hide

# We can also solve this problem with the native [Boris pusher](@ref demo_boris).
# Note that the Boris pusher requires a additional parameters: a fixed timestep, and an output save interval.

dt = 0.1
savestepinterval = 1

## Solve for the trajectories

prob = TraceProblem(stateinit, tspan, param; prob_func)
trajs = TestParticle.solve(prob; dt, trajectories, savestepinterval)

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Electron trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data,
)

for i in eachindex(trajs)
   lines!(ax, trajs[i]; idxs=(1,2,3), label="$i")
   ##TODO: wait for https://github.com/MakieOrg/Makie.jl/issues/3623 to be fixed!
   ax.scene.plots[9+2*i-1].color = Makie.wong_colors()[i]
end

f = DisplayAs.PNG(f) #hide

# Note that by default linear interpolation is applied when plotting the trajectories from the Boris method.