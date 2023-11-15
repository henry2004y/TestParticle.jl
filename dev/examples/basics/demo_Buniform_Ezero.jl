import DisplayAs #hide
using TestParticle
using TestParticle: get_gc
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png")

uniform_B(x) = SA[0.0, 0.0, 1e-8]
uniform_E(x) = SA[0.0, 0.0, 0.0]

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 18)

param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Tsit5())

f = plot(sol)
f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
