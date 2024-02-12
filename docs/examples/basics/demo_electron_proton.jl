# ---
# title: Proton and electron in a static EM field
# id: demo_proton_electron
# date: 2023-04-20
# author: "[Hongyang Zhou](https://github.com/henry2004y); [Tiancheng Liu](https://github.com/TCLiuu)"
# julia: 1.9.0
# description: Simple particle trajectories under uniform B and zero E
# ---

# This example demonstrates tracing a single proton and electron motion under a uniform B field in real physical parameters. The E field is assumed to be zero such that there is no particle acceleration.
# Due to the fact that ``m_p / m_e \doteq 1836``, the proton gyro-radius is 1800 times larger than the electron, if they start with the same velocity as in this case. In more common cases we would compare electrons and protons with the same energy, and their gyro-radii differ by a factor of ``\sqrt{m_p/m_e} \sim 40``.

import DisplayAs #hide
using TestParticle
using Meshes
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png") #hide

### Initialize grid and field

x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
z = range(-10, 10, length=25)
B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

B[3,:,:,:] .= 1e-11
E[3,:,:,:] .= 5e-13

Δx = x[2] - x[1]
Δy = y[2] - y[1]
Δz = z[2] - z[1]

grid = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
   (x[1], y[1], z[1]),
   (Δx, Δy, Δz))

### Initialize particles

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param_electron = prepare(grid, E, B, species=Electron)
tspan_electron = (0.0, 15.0)

param_proton = prepare(grid, E, B, species=Proton)
tspan_proton = (0.0, 10.0)

### Solve for the trajectories

prob_e = ODEProblem(trace!, stateinit, tspan_electron, param_electron)
prob_p = ODEProblem(trace!, stateinit, tspan_proton, param_proton)

sol_e = solve(prob_e, Vern9())
sol_p = solve(prob_p, Vern9())

### Visualization

f = Figure()
Axis3(f[1,1], aspect = :data)
plot!(sol_e, color=:tomato, label="electron")
plot!(sol_p, color=:deepskyblue3, label="proton")
axislegend()

f = DisplayAs.PNG(f) #hide