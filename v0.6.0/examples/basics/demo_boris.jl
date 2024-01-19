import DisplayAs #hide
using TestParticle
using StaticArrays
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png")

uniform_B(x) = SA[0.0, 0.0, 0.01]
uniform_E(x) = SA[0.0, 0.0, 0.0]

x0 = [0.0, 0.0, 0.0]
v0 = [0.0, 1e5, 0.0]
stateinit = [x0..., v0...]
tspan = (0.0, 3e-6)
dt = 3e-11
param = prepare(uniform_E, uniform_B, species=Electron)
paramBoris = BorisMethod(param)
prob = TraceProblem(stateinit, tspan, dt, paramBoris)

traj = trace_trajectory(prob; savestepinterval=10);
@time traj = trace_trajectory(prob; savestepinterval=10);

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=10*dt);
@time sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=10*dt);

sol2 = solve(prob, Tsit5());
@time sol2 = solve(prob, Tsit5());

f = Figure(size=(700, 600))
ax = Axis(f[1, 1], aspect=1)
@views lines!(ax, traj[1,:], traj[2,:], label="Boris")
lines!(ax, sol1, linestyle=:dashdot, label="Tsit5 fixed")
lines!(ax, sol2, linestyle=:dot, label="Tsit5 adaptive")
axislegend(position=:lt, framevisible=false)
f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
