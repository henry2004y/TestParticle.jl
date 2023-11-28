import DisplayAs #hide
using TestParticle
using TestParticle: get_gc
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png")

B(x) = SA[0.0, 1e-8, 0.0]
E(x) = SA[0.0, 0.0, 0.0]

# Earth's gravity
F(x) = SA[0.0, 0.0, -TestParticle.mᵢ*9.8]

# initial static particle
x0 = [1.0, 0, 0]
v0 = [0.0, 0.0, 0.0]
stateinit = [x0..., v0...]
tspan = (0, 1.0)

param = prepare(E, B, F, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())
# drift in x-direction + free fall in z-direction
f = orbit(sol, vars=[(1,3),], interactive=false)

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
