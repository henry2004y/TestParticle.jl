# ---
# title: Helix motion
# id: demo_uniformB_zeroE
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Simple proton trajectory under uniform B and zero E
# ---

# This example demonstrates a single proton motion under a uniform B field. The E field is assumed to be zero such that there is no particle acceleration.

import DisplayAs #hide
using TestParticle
using TestParticle: get_gc
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png") #hide

uniform_B(x) = SA[0.0, 0.0, 1e-8]
uniform_E(x) = SA[0.0, 0.0, 0.0]

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 18)

param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

f = plot(sol)
f = DisplayAs.PNG(f) #hide