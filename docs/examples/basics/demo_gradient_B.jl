# ---
# title: Grad-B drift
# id: demo_gradB
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Simple magnetic field gradient drift demonstration using Makie
# ---

# This example demonstrates a single proton motion under a non-uniform B field with gradient ∇B ⊥ B.
# The orbit of guiding center includes some high order terms, it is different from the
# formula of magnetic field gradient drift of some textbooks which just preserves the first
# order term.
# More theoretical details can be found in Introduction to Plasma Physics and Controlled Fusion by F. F. Chen, and Fundamentals of Plasma Physics by Paul Bellan.

import DisplayAs # hide

using TestParticle
using TestParticle: get_gc
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using ForwardDiff: gradient
using CairoMakie
CairoMakie.activate!(type = "png")

function grad_B(x)
    return SA[0, 0, 1e-8+1e-9 *x[2]]
end

function uniform_E(x)
    return SA[1e-9, 0, 0]
end

abs_B(x) = norm(grad_B(x))

## trace the orbit of the guiding center
function trace_gc!(dx, x, p, t)
    q, m, E, B, sol = p
    xu = sol(t)
    gradient_B = gradient(abs_B, x)
    Bv = B(x)
    b = normalize(Bv)
    v_par = (xu[4:6]⋅b).*b
    v_perp = xu[4:6] - v_par
    dx[1:3] = m*norm(v_perp)^2*(Bv×gradient_B)/(2*q*norm(Bv)^3) + (E(x)×Bv)/norm(Bv)^2+v_par
end

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 20)
param = prepare(uniform_E, grad_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Tsit5(); save_idxs=[1,2,3,4,5,6])

gc = get_gc(param)
gc_x0 = [gc_i(stateinit) for gc_i in gc]
prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Tsit5(); save_idxs=[1,2,3])

gc_analytic = Tuple(xu -> getindex(sol_gc(xu[7]), i) for i = 1:3)
## numeric result and analytic result
f = orbit(sol, vars=[(1, 2, 3), gc, gc_analytic])

f = DisplayAs.PNG(f) # hide