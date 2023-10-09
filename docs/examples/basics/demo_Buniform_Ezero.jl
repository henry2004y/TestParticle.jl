# ---
# title: Helix motion
# id: demo_uniformB_zeroE
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Simple proton trajectory under uniform B and zero E
# ---

# This example demonstrates a single proton motion under a uniform B field. The E field is assumed to be zero such that there is no particle acceleration.

#using JSServe: Page # hide
#Page(exportable=true, offline=true) # hide
#import DisplayAs # hide

using TestParticle
using TestParticle: get_gc
using TestParticleMakie
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
#using WGLMakie

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
#f = DisplayAs.PNG(f)