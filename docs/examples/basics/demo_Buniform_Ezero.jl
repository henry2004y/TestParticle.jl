# ---
# title: Helix motion
# id: demo_uniformB_zeroE
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Simple proton trajectory under uniform B and zero E
# ---

# This example demonstrates a single proton motion under a uniform B field.
# The E field is assumed to be zero such that there is no particle acceleration.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

uniform_B(x) = SA[0.0, 0.0, 1e-8]
uniform_E(x) = SA[0.0, 0.0, 0.0]

## Initial condition
stateinit = let x0 = [1.0, 0, 0], v0 = [0.0, 1.0, 0.1]
   [x0..., v0...]
end
## Time span
tspan = (0, 18)

param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

## Visualization
f = Figure(fontsize=18)
ax = Axis3(f[1, 1],
   title = "Helix Trajectory",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data,
)

plot!(ax, sol, idxs=(1, 2, 3))

f = DisplayAs.PNG(f) #hide