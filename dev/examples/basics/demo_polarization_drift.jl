import DisplayAs #hide
using TestParticle
using TestParticle: get_gc
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png")

function uniform_B(x)
    return SA[0, 0, 1e-8]
end

function time_varying_E(x, t)
    # return SA[0, 1e-9*cos(0.1*t), 0]
    return SA[0, 1e-9*0.1*t, 0]
end

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 100)
param = prepare(time_varying_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Tsit5(); save_idxs=[1,2,3,4,5,6])
gc = get_gc(param)
v_perp(xu) = hypot(xu[4], xu[5])
gc_y = gc[2]
# polarization drift
f = monitor(sol, vars=[v_perp, 2, gc_y])

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl