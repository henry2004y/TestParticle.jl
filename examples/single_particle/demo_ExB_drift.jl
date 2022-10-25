# E×B drift
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

# trace the orbit of the guiding center
function trace_gc!(dx, x, p, t)
    _, _, E, B, sol = p
    xu = sol(t)
    Bv = B(x)
    b = normalize(Bv)
    v_para = (xu[4:6]⋅b).*b
    B2 = sum(Bv.^2)
    dx[1:3] = (E(x)×Bv)/B2 + v_para
end

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 20)
# E×B drift
param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Tsit5(); save_idxs=[1,2,3,4,5,6])

gc = get_gc(param)
gc_x0 = [gc_i(stateinit) for gc_i in gc]
prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Tsit5(); save_idxs=[1,2,3])

gc_analytic = Tuple(xu -> getindex(sol_gc(xu[7]), i) for i = 1:3)
# numeric result and analytic result
orbit(sol, vars=[(1, 2, 3), gc, gc_analytic])