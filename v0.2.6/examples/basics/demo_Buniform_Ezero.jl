using JSServe: Page # hide
Page(exportable=true, offline=true) # hide

using TestParticle
using TestParticle: get_gc
using TestParticleMakie
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using WGLMakie

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

plot(sol)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
