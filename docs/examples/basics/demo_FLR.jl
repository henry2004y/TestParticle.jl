# ---
# title: Finite-Larmor-Radius effect
# id: demo_FLR
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu), [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.0
# description: Finite Larmor radius effect demonstration using Makie
# ---

# The general FLR effect refers to the correction terms introduced when considering the field difference at the particle location and the guiding center location.
# More theoretical details can be found in [Non-uniform E Field](https://henry2004y.github.io/KeyNotes/contents/single.html#sec-nonuniform_E).

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra: ×, ⋅, norm, normalize
using Tensors: laplace
import Tensors: Vec as Vec3
## using SpecialFunctions
using CairoMakie
CairoMakie.activate!(type = "png") #hide

function uniform_B(x)
    return SA[0, 0, 1e-8]
end

function nonuniform_E(x)
    return SA[1e-9*cos(0.3*x[1]), 0, 0]
end

## trace the orbit of the guiding center
function trace_gc!(dx, x, p, t)
    q2m, E, B, sol = p
    xu = sol(t)
    xp = @view xu[1:3]
    Bv = B(xp)
    b = normalize(Bv)
    v_par = (xu[4:6] ⋅ b) .* b # (v ⋅ b)b
    v_perp = xu[4:6] - v_par
    r4 = (norm(v_perp) / q2m / norm(Bv))^2 / 4
    EB(x) = (E(x) × B(x)) / norm(B(x))^2
    ## dx[1:3] = EB(xp) + v_par
    dx[1:3] = EB(x) + r4*laplace.(EB, Vec3(x...)) + v_par

    ## more accurate
    ## dx[1:3] = besselj0(0.3*norm(v_perp)/q2m/norm(Bv))*EB(x) + v_par
end

## Initial condition
stateinit = let x0 = [1.0, 0, 0], v0 = [0.0, 1.0, 0.1]
    [x0..., v0...]
end
## Time span
tspan = (0, 20)
param = prepare(nonuniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

gc = get_gc(param)
gc_x0 = gc(stateinit)
prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Tsit5(); save_idxs=[1,2,3])

## numeric result and analytic result
f = Figure(fontsize=18)
ax = Axis3(f[1, 1],
   title = "Finite Larmor Radius Effect",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data,
   azimuth = 0.3π,
)

gc_plot(x, y, z, vx, vy, vz) = (gc(SA[x, y, z, vx, vy, vz])...,)

lines!(ax, sol, idxs=(1, 2, 3))
lines!(ax, sol, idxs=(gc_plot, 1, 2, 3, 4, 5, 6))
lines!(ax, sol_gc, idxs=(1, 2, 3))

for i in 1:3
    ##TODO: wait for https://github.com/MakieOrg/Makie.jl/issues/3623 to be fixed!
    ax.scene.plots[9+i].color = Makie.wong_colors()[i]
end

f = DisplayAs.PNG(f) #hide