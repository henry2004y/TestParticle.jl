# ---
# title: E×B drift
# id: demo_ExB
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Simple ExB drift demonstration
# ---

# This example demonstrates a single proton motion under uniform E and B fields. The electric field is parallel to the magnetic field in the z-direction, so the motion consists of a cyclotron gyration and an acceleration along z. On top of that, particles also exhibit a ExB drift in the direction perpendicular to both E and B field.
# More theoretical details can be found in [ExB Drift](https://henry2004y.github.io/KeyNotes/contents/single.html#finite-e) and Computational Plasma Physics, Toshi Tajima.

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

function uniform_E(x)
    return SA[1e-9, 0, 0]
end

## Trace the orbit of the guiding center
function trace_gc!(dx, x, p, t)
    _, E, B, sol = p
    xu = sol(t)
    Bv = B(x)
    b = normalize(Bv)
    v_par = (xu[4:6] ⋅ b) .* b
    B2 = sum(Bv.^2)
    dx[1:3] = (E(x) × Bv)/B2 + v_par
end

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 20)
## E×B drift
param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())
## Functions for obtaining the guiding center from actual trajectory
gc = get_gc(param)
gc_x0 = [gc_i(stateinit) for gc_i in gc]
prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Vern9(); save_idxs=[1,2,3])
## Functions for obtaining the guiding center from analytic drifts
gc_analytic = Tuple(xu -> getindex(sol_gc(xu[7]), i) for i = 1:3)
## Numeric and analytic results
f = orbit(sol, vars=[(1, 2, 3), gc, gc_analytic])

f = DisplayAs.PNG(f) #hide

# Note that in this simple ExB drift case, the analytic and numeric guiding centers have the same motion.