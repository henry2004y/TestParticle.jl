using JSServe: Page # hide
Page(exportable=true, offline=true) # hide

using TestParticle
using TestParticle: get_gc
using TestParticleMakie
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using WGLMakie
using Tensors: laplace
import Tensors: Vec as Vec3
# using SpecialFunctions

function uniform_B(x)
    return SA[0, 0, 1e-8]
end

function nonuniform_E(x)
    return SA[1e-9*cos(0.3*x[1]), 0, 0]
end

# trace the orbit of the guiding center
function trace_gc!(dx, x, p, t)
    q, m, E, B, sol = p
    xu = sol(t)
    xp = @view xu[1:3]
    Bv = B(xp)
    b = normalize(Bv)
    v_par = (xu[4:6]⋅b).*b  # (v⋅b)b
    v_perp = xu[4:6] - v_par
    r4 = (m*norm(v_perp)/q/norm(Bv))^2/4
    EB(x) = (E(x)×B(x))/norm(B(x))^2
    # dx[1:3] = EB(xp) + v_par
    dx[1:3] = EB(x) + r4*laplace.(EB, Vec3(x...)) + v_par

    # more accurate
    # dx[1:3] = besselj0(0.3*m*norm(v_perp)/q/norm(Bv))*EB(x) + v_par
end

x0 = [1.0, 0, 0]
v0 = [0.0, 1.0, 0.1]
stateinit = [x0..., v0...]
tspan = (0, 20)
param = prepare(nonuniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Tsit5(); save_idxs=[1,2,3,4,5,6])

gc = get_gc(param)
gc_x0 = [gc_i(stateinit) for gc_i in gc]
prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Tsit5(); save_idxs=[1,2,3])

gc_analytic = Tuple(xu -> getindex(sol_gc(xu[7]), i) for i = 1:3)
# numeric result and analytic result
orbit(sol, vars=[(1, 2, 3), gc, gc_analytic])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
