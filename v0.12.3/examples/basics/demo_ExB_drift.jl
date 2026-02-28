import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra: ⋅, ×, normalize
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# Analytic EM fields
uniform_B(x) = SA[0, 0, 1e-8]
uniform_E(x) = SA[1e-9, 0, 0]

# Trace the orbit of the guiding center
function trace_gc_ExB!(dx, x, p, t)
   _, E, B, sol = p
   xu = sol(t)
   Bv = B(x)
   b = normalize(Bv)
   v_par = @views (xu[4:6] ⋅ b) .* b
   B2 = sum(Bv.^2)
   dx[1:3] = (E(x) × Bv) / B2 + v_par
end
# Initial condition
stateinit = let x0 = [1.0, 0.0, 0.0], v0 = [0.0, 1.0, 0.1]
   [x0..., v0...]
end
# Time span
tspan = (0, 20)

# Trace particle
param = prepare(uniform_E, uniform_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

# Functions for obtaining the guiding center from actual trajectory
gc = get_gc(param)
gc_x0 = gc(stateinit)
prob_gc = ODEProblem(trace_gc_ExB!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Vern9(); save_idxs=[1,2,3]);

# Numeric and analytic results
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

lines!(ax, sol, idxs=(1, 2, 3), color=Makie.wong_colors()[1])
lines!(ax, sol, idxs=(gc_plot, 1, 2, 3, 4, 5, 6), color=Makie.wong_colors()[2])
lines!(ax, sol_gc, idxs=(1, 2, 3), linestyle=:dash, color=Makie.wong_colors()[3])

f = DisplayAs.PNG(f) #hide

stateinit_gc, param_gc = prepare_gc(stateinit, uniform_E, uniform_B, species=Proton, removeExB=false)
prob_gc1st = ODEProblem(trace_gc_1st!, stateinit_gc, tspan, param_gc)
sol_gc1st = solve(prob_gc1st, Vern9())

f = Figure(fontsize=18)
ax = Axis3(f[1, 1],
   title = "ExB Drift",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data,
   azimuth = 0.3π,
)

lines!(ax, sol, idxs=(1, 2, 3), color=Makie.wong_colors()[1])
lines!(ax, sol_gc1st, idxs=(1, 2, 3), linestyle=:dash, color=Makie.wong_colors()[3])

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
