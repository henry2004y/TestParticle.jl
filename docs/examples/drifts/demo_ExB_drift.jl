# # E×B Drift
#
# This example demonstrates a single proton motion under uniform E and B fields.
# The electric field is parallel to the magnetic field in the z-direction, so the motion consists of a cyclotron gyration and an acceleration along z.
# On top of that, particles also exhibit an ExB drift in the direction perpendicular to both E and B field.
# More theoretical details can be found in [ExB Drift](https://henry2004y.github.io/KeyNotes/contents/single.html#finite-e).

import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using LinearAlgebra: ⋅, ×, normalize
using TestParticle: get_EField, get_BField
using CairoMakie
CairoMakie.activate!(type = "png") #hide

## Analytic EM fields
uniform_B(x) = SA[0, 0, 1e-8]
uniform_E(x) = SA[1e-9, 0, 0]

## Initial condition
stateinit = let x0 = [1.0, 0.0, 0.0], v0 = [0.0, 1.0, 0.1]
   [x0..., v0...]
end
## Time span
tspan = (0, 20)

## Trace particle
param = prepare(uniform_E, uniform_B, species = Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

## Functions for obtaining the guiding center from actual trajectory
gc = param |> get_gc_func
gc_x0 = gc(stateinit) |> Vector
prob_gc = ODEProblem(trace_gc_exb!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Vern9(); save_idxs = [1, 2, 3]);

## Numeric and analytic results
f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "ExB Drift",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data,
   azimuth = 0.3π
)

gc_plot(x, y, z, vx, vy, vz) = (gc(SA[x, y, z, vx, vy, vz])...,)

lines!(ax, sol, idxs = (1, 2, 3), color = Makie.wong_colors()[1])
lines!(ax, sol, idxs = (gc_plot, 1, 2, 3, 4, 5, 6), color = Makie.wong_colors()[2])
lines!(ax, sol_gc, idxs = (1, 2, 3), linestyle = :dash, color = Makie.wong_colors()[3])

f = DisplayAs.PNG(f) #hide

# Note that in this simple ExB drift case, the analytic and numeric guiding centers overlaps. Also note that `trace_gc_exb!` here depends on the velocity at time `t` from the particle trajectory, which is not exactly the guiding center velocity.
# A first-order GC approximation tracker would be the following:

stateinit_gc,
param_gc = prepare_gc(stateinit, uniform_E, uniform_B, species = Proton, removeExB = false)
prob_gc1st = ODEProblem(trace_gc_1st!, stateinit_gc, tspan, param_gc)
sol_gc1st = solve(prob_gc1st, Vern9())

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "ExB Drift",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data,
   azimuth = 0.3π
)

lines!(ax, sol, idxs = (1, 2, 3), color = Makie.wong_colors()[1])
lines!(ax, sol_gc1st, idxs = (1, 2, 3), linestyle = :dash, color = Makie.wong_colors()[3])

f = DisplayAs.PNG(f) #hide
