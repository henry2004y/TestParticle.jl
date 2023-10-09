import DisplayAs # hide

using TestParticle
using TestParticle: get_gc
using TestParticleMakie
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png")

function uniform_B(x)
    return SA[0.0, 0.0, 1e-8]
end

function uniform_E(x)
    return SA[0.0, 0.0, 0.0]
end

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 18)

param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Tsit5())

f = plot(sol)
f = DisplayAs.PNG(f) # hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
