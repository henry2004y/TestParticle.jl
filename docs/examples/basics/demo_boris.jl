# ---
# title: Boris method
# id: demo_boris
# date: 2023-11-13
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.0
# description: Simple electron trajectory under uniform B and zero E
# ---

# This example demonstrates a single electron motion under a uniform B field. The E field is assumed to be zero such that there is no particle acceleration.
# We use the [Boris method](https://apps.dtic.mil/sti/citations/ADA023511) for phase space conservation under a fixed time step.
# This is compared against other ODE general algorithms for performance and accuracy.

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

# Now let's compare against other ODE solvers:

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=10*dt);
@time sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=10*dt);

sol2 = solve(prob, Tsit5());
@time sol2 = solve(prob, Tsit5());

# The phase space conservation can be visually checked by plotting the trajectory:

f = Figure(size=(700, 600))
ax = Axis(f[1, 1], aspect=1)
@views lines!(ax, traj[1][1,:], traj[1][2,:], label="Boris")
lines!(ax, sol1, linestyle=:dashdot, label="Tsit5 fixed")
lines!(ax, sol2, linestyle=:dot, label="Tsit5 adaptive")
axislegend(position=:lt, framevisible=false)
f = DisplayAs.PNG(f) #hide

# We can see that the Boris method is both accurate and fast; Fixed time step `Tsit5()` is ok, but adaptive `Tsit5()` is pretty bad for long time evolutions.
# When calling OrdinaryDiffEq.jl, we recommend using `Vern9()` as a starting point instead of `Tsit5()`, especially combined with adaptive timestepping.