import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

"""
Set initial state for EnsembleProblem.
"""
function prob_func(prob, i, repeat)
   prob = @views remake(prob, u0 = [prob.u0[1:3]..., i/3, 0.0, 0.0])
end

# Initialization

B(x) = SA[0, 0, 1e-11]
E(x) = SA[0, 0, 1e-13]

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param = prepare(E, B, species = Electron)
tspan = (0.0, 10.0)

trajectories = 3

# Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy = false)
sols = solve(ensemble_prob, Vern7(), EnsembleThreads(); trajectories)

# Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Electron trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data
)

for i in eachindex(sols)
   lines!(ax, sols[i], idxs = (1, 2, 3), label = "$i", color = Makie.wong_colors()[i])
end

f = DisplayAs.PNG(f) #hide

dt = 0.1
savestepinterval = 1

# Solve for the trajectories

prob = TraceProblem(stateinit, tspan, param; prob_func)
trajs = TestParticle.solve(prob; dt, trajectories, savestepinterval)

# Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Electron trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data
)

for i in eachindex(trajs)
   lines!(ax, trajs[i]; idxs = (1, 2, 3), label = "$i", color = Makie.wong_colors()[i])
end

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
