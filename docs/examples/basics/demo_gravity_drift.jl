# ---
# title: ExF drift
# id: demo_uniformB_gravity
# date: 2023-04-19
# author: "[Hongyang Zhou](https://github.com/henry2004y); [Tiancheng Liu](https://github.com/TCLiuu)"
# julia: 1.9.0
# description: Simple proton trajectory under uniform B and gravity
# ---

# This example demonstrates a single proton motion under uniform B and gravity fields.

import DisplayAs #hide
using TestParticle
using TestParticle: get_gc
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

B(x) = SA[0.0, 1e-8, 0.0]
E(x) = SA[0.0, 0.0, 0.0]

## Earth's gravity
F(x) = SA[0.0, 0.0, -TestParticle.mᵢ*9.8]

## initial static particle
x0 = [1.0, 0, 0]
v0 = [0.0, 0.0, 0.0]
stateinit = [x0..., v0...]
tspan = (0, 1.0)

param = prepare(E, B, F, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())
## drift in x-direction + free fall in z-direction
f = lines(sol, idxs=(3,1);
   figure = (; size = (800, 400), fontsize=18),
   axis = (; title="ExF Drift", xlabel="Z [m]", ylabel="X [m]", aspect = DataAspect())
)

f = DisplayAs.PNG(f) #hide