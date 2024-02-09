# ---
# title: Boris method
# id: demo_boris
# date: 2024-01-23
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

## Reference parameters
const tperiod = 2π / (abs(param[1]) * hypot(param[3]([0.0,0.0,0.0], 0.0)...))
const rL = hypot(v0...) / (abs(param[1]) * Bmag)
const invrL = 1 / rL;

# We first trace the particle for one period with a discrete time step of a quarter period.

tspan = (0.0, tperiod)
dt = tperiod / 4

prob = TraceProblem(stateinit, tspan, param)

sol_boris = TestParticle.solve(prob; dt, savestepinterval=1);

# Let's compare against the default ODE solver `Tsit5` from DifferentialEquations.jl, in both fixed time step mode and adaptive mode:

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);
sol2 = solve(prob, Tsit5());

## Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

# It is clear that the Boris method comes with larger phase errors (``\mathcal{O}(\Delta t)``) compared with Tsit5.
# The phase error gets smaller using a smaller dt:

dt = tperiod / 8

prob = TraceProblem(stateinit, tspan, param)

sol_boris = TestParticle.solve(prob; dt, savestepinterval=1);

prob = ODEProblem(trace!, stateinit, tspan, param)
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);

## Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

# The Boris pusher shines when we do long time tracing, which is fast and conserves energy:

tspan = (0.0, 200*tperiod)
dt = tperiod / 12

prob_boris = TraceProblem(stateinit, tspan, param)
prob = ODEProblem(trace!, stateinit, tspan, param)

sol_boris = TestParticle.solve(prob_boris; dt, savestepinterval=10);
sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);
sol2 = solve(prob, Tsit5());

## Visualization
f = plot_trajectory(sol_boris, sol1, sol2)
f = DisplayAs.PNG(f) #hide

# Fixed time step `Tsit5()` is ok, but adaptive `Tsit5()` is pretty bad for long time evolutions. The change in radius indicates change in energy, which is sometimes known as numerical heating.
#
# Another aspect to compare is performance:

@time sol_boris = TestParticle.solve(prob_boris; dt, savestepinterval=10);
@time sol1 = solve(prob, Tsit5(); adaptive=false, dt, dense=false, saveat=dt);
@time sol2 = solve(prob, Tsit5());

# We can extract the solution `(x, y, z, vx, vy, vz)` at any given time by performing a linear interpolation:

t = tspan[2] / 2
sol_boris[1](t)

# The Boris method is faster and consumes less memory. However, in practice, it is pretty hard to find an optimal algorithm.
# When calling OrdinaryDiffEq.jl, we recommend using `Vern9()` as a starting point instead of `Tsit5()`, especially combined with adaptive timestepping.