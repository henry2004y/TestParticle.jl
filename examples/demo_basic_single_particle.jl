# Some basic single particle cases.
# More theoretical details can be found in F.F.Chen's Introduction to Plasma Physics and Controlled Fusion.

using TestParticle
using TestParticle: get_gc
using TestParticleMakie
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using GLMakie

function uniform_B(x)
    return SA[0, 0, 1e-8]
end

function uniform_E(x)
    return SA[1e-9, 0, 0]
end

function zero_E(x)
    return SA[0, 0, 0]
end

function grad_B(x)
    return SA[0, 0, 1e-8+1e-9 *x[2]]
end

function curved_B(x)
    # satisify ∇⋅B=0
    θ = atan(x[3]/(x[1]+3))
    return SA[-1e-8*sin(θ), 0, 1e-8*cos(θ)]
end

function nonuniform_E(x)
    return SA[1e-9*cos(2*π*x[1]), 0, 0]
end

function time_varying_E(x, t)
    # return SA[0, 1e-9*cos(0.1*t), 0]
    return SA[0, 1e-9*0.1*t, 0]
end

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 20)

# E×B drift
param = prepare(uniform_E, uniform_B, species=Proton)
prob_EB = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob_EB, Tsit5(); save_idxs=[1,2,3,4,5,6])
gc = get_gc(param)
orbit(sol, vars=[(1, 2, 3), gc])

# magnetic field gradient drift
param = prepare(zero_E, grad_B, species=Proton)
prob_grad = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob_grad, Vern8(); save_idxs=[1,2,3,4,5,6])
gc = get_gc(param)
# The orbit of guiding center includes some high order terms, it is different from the formula of
# magnetic field gradient drift of some textbooks which just preserves the first order term.
orbit(sol, vars=[(1, 2, 3), gc])

# magnetic field curvature drift
param = prepare(zero_E, curved_B, species=Proton)
prob_curve = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob_curve, Vern8(); save_idxs=[1,2,3,4,5,6])
gc = get_gc(param)
orbit(sol, vars=[(1, 2, 3), gc])

# finite Larmor radius effect
param = prepare(nonuniform_E, uniform_B, species=Proton)
prob_FLR = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob_FLR, Vern8(); save_idxs=[1,2,3,4,5,6])
gc = get_gc(param)
orbit(sol, vars=[(1, 2, 3), gc])

# polarization drift
param = prepare(time_varying_E, uniform_B, species=Proton)
tspan = (0, 200)
prob_pd = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob_pd, Vern8(); save_idxs=[1,2,3,4,5,6])
gc = get_gc(param)
v_perp(xu) = hypot(xu[4], xu[5])
gc_y = gc[2]
monitor(sol, vars=[v_perp, 2, gc_y])