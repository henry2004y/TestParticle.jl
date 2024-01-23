import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png")

"Set initial state for EnsembleProblem."
function prob_func(prob, i, repeat)
   # If u0 is immutable (e.g. StaticArrays), use remake; otherwise we can directly modify u0
   ##remake(prob, u0=[i/3, 0.0, 0.0])
   prob.u0[4] = i / 3

   prob
end

# Initialization

B(x) = SA[0, 0, 1e-11]
E(x) = SA[0, 0, 1e-13]

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param = prepare(E, B, species=Electron)
tspan = (0.0, 10.0)

trajectories = 3

# Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
sols = solve(ensemble_prob, Tsit5(), EnsembleThreads(); trajectories)

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

f = DisplayAs.PNG(f) #hide

dt = 0.1
savestepinterval = 1

# Solve for the trajectories

paramBoris = BorisMethod(param)
prob = TraceProblem(stateinit, tspan, dt, paramBoris; prob_func)
trajs = trace_trajectory(prob; trajectories, savestepinterval)

# Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Electron trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data,
)

for i in eachindex(trajs)
   @views lines!(ax, trajs[i][1,:], trajs[i][2,:], trajs[i][3,:], label="$i")
end

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
