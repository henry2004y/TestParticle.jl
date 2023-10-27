import DisplayAs # hide

using TestParticle
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png")

"Initial state perturbation for EnsembleProblem."
function prob_func(prob, i, repeat)
   remake(prob, u0=rand()*prob.u0)
end

# Initialization

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

# Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sols = solve(ensemble_prob, Tsit5(), EnsembleThreads();
   trajectories=trajectories, save_idxs=[1,2,3])

# Visualization

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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
