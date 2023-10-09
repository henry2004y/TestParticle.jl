import DisplayAs # hide

using TestParticle
using TestParticle: get_gc
using TestParticleMakie
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png")

function B(x)
    return SA[0.0, 1e-8, 0.0]
end

function E(x)
    return SA[0.0, 0.0, 0.0]
end
# Earth's gravity
function F(x)
    return SA[0.0, 0.0, -TestParticle.mᵢ*9.8]
end
# initial static particle
x0 = [1.0, 0, 0]
v0 = [0.0, 0.0, 0.0]
stateinit = [x0..., v0...]
tspan = (0, 1.0)

param = prepare(E, B, F, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Tsit5(); save_idxs=[1,2,3])
# drift in x-direction + free fall in z-direction
f = plot(sol)

f = DisplayAs.PNG(f) # hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
