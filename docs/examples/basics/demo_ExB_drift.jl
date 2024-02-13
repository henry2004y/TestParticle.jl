# ---
# title: E×B drift
# id: demo_ExB
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.0
# description: Simple ExB drift demonstration
# ---

# This example demonstrates a single proton motion under uniform E and B fields.
# The electric field is parallel to the magnetic field in the z-direction, so the motion consists of a cyclotron gyration and an acceleration along z.
# On top of that, particles also exhibit an ExB drift in the direction perpendicular to both E and B field.
# Note that in this simple ExB drift case, the analytic and numeric guiding centers overlaps.
# More theoretical details can be found in [ExB Drift](https://henry2004y.github.io/KeyNotes/contents/single.html#finite-e).

import DisplayAs #hide
using TestParticle
using TestParticle: get_gc
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra: ⋅, ×, normalize
using CairoMakie
CairoMakie.activate!(type = "png") #hide

## Analytic EM fields
uniform_B(x) = SA[0, 0, 1e-8]
uniform_E(x) = SA[1e-9, 0, 0]

## Trace the orbit of the guiding center
function trace_gc!(dx, x, p, t)
    _, E, B, sol = p
    xu = sol(t)
    Bv = B(x)
    b = normalize(Bv)
    v_par = @views (xu[4:6] ⋅ b) .* b
    B2 = sum(Bv.^2)
    dx[1:3] = (E(x) × Bv) / B2 + v_par
end
## Initial condition
stateinit = let x0 = [1.0, 0.0, 0.0], v0 = [0.0, 1.0, 0.1]
    [x0..., v0...] 
end
## Time span
tspan = (0, 20)

## Trace particle
param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

## Functions for obtaining the guiding center from actual trajectory
gc = get_gc(param)
gc_x0 = gc(stateinit)
prob_gc = ODEProblem(trace_gc!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Vern9(); save_idxs=[1,2,3]);

## Numeric and analytic results
f = Figure(fontsize=18)
ax = Axis3(f[1, 1],
   title = "ExB Drift",
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
    ax.scene.plots[9+2*i-1].color = Makie.wong_colors()[i]
end

f = DisplayAs.PNG(f) #hide