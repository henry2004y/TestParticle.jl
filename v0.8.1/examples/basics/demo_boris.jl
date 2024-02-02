import DisplayAs #hide
using TestParticle
using StaticArrays
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png")

function plot_trajectory(sol_boris, sol1, sol2)
   f = Figure(size=(700, 600))
   ax = Axis(f[1, 1], aspect=1, limits = (-3, 1, -2, 2))
   @views lines!(ax, sol_boris[1].u[1,:]./rL, sol_boris[1].u[2,:]./rL,
      linewidth=2, label="Boris")
   l1 = lines!(ax, sol1, linestyle=:dashdot, label="Tsit5 fixed", linewidth=2)
   l2 = lines!(ax, sol2, linestyle=:dot, label="Tsit5 adaptive", linewidth=2)
   scale!(l1, invrL, invrL, invrL)
   scale!(l2, invrL, invrL, invrL)
   axislegend(position=:lt, framevisible=false)

   f
end

const Bmag = 0.01
uniform_B(x) = SA[0.0, 0.0, Bmag]
uniform_E(x) = SA[0.0, 0.0, 0.0]

x0 = [0.0, 0.0, 0.0]
v0 = [0.0, 1e5, 0.0]
stateinit = [x0..., v0...]

param = prepare(uniform_E, uniform_B, species=Electron)

# Reference parameters
const tperiod = 2π / (abs(param[1]) * hypot(param[3]([0.0,0.0,0.0], 0.0)...))
const rL = hypot(v0...) / (abs(param[1]) * Bmag)
const invrL = 1 / rL;

tspan = (0.0, tperiod)
dt = tperiod / 4

prob = TraceProblem(stateinit, tspan, dt, param)

sol_boris = TestParticle.solve(prob; savestepinterval=1);

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);
sol2 = solve(prob, Tsit5());

# Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

dt = tperiod / 8

prob = TraceProblem(stateinit, tspan, dt, param)

sol_boris = TestParticle.solve(prob; savestepinterval=1);

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);

# Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

tspan = (0.0, 200*tperiod)
dt = tperiod / 12

prob_boris = TraceProblem(stateinit, tspan, dt, param)
prob = ODEProblem(trace!, stateinit, tspan, param)

sol_boris = TestParticle.solve(prob_boris; savestepinterval=10);
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);
sol2 = solve(prob, Tsit5());

# Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

@time sol_boris = TestParticle.solve(prob_boris; savestepinterval=10);
@time sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);
@time sol2 = solve(prob, Tsit5());

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
