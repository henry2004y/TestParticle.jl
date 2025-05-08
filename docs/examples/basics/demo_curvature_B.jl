# ---
# title: Curl-B drift
# id: demo_curlB
# date: 2023-04-19
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Simple magnetic field curvature drift demonstration
# ---

# This example demonstrates a single proton motion under a vacuum non-uniform B field with gradient and curvature.
# The analytic calculation includes the grad-B drift, the curvature drift, the ExB drift and parallel velocity.
# More theoretical details can be found in [Curvature Drift](https://henry2004y.github.io/KeyNotes/contents/single.html#curved-b-curvature-drift) and Computational Plasma Physics by Toshi Tajima.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using LinearAlgebra: normalize, norm, ×, ⋅
using ForwardDiff: gradient, jacobian
using CairoMakie
CairoMakie.activate!(type = "png") #hide

function curved_B(x)
    ## satisify ∇⋅B=0
    ## B_θ = 1/r => ∂B_θ/∂θ = 0
    θ = atan(x[3] / (x[1] + 3))
    r = sqrt((x[1] + 3)^2 + x[3]^2)
    return SA[-1e-7*sin(θ)/r, 0, 1e-7*cos(θ)/r]
end

function zero_E(x)
    return SA[0, 0, 0]
end

abs_B(x) = norm(curved_B(x))  # |B|

## Initial conditions
stateinit = let x0 = [1.0, 0.0, 0.0], v0 = [0.0, 1.0, 0.1]
    [x0..., v0...]
end
## Time span
tspan = (0, 40)
## Trace particle
param = prepare(zero_E, curved_B, species=Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())
## Functions for obtaining the guiding center from actual trajectory
gc = param |> get_gc_func
gc_x0 = gc(stateinit) |> Vector
prob_gc = ODEProblem(trace_gc_drifts!, gc_x0, tspan, (param..., sol))
sol_gc = solve(prob_gc, Vern7(); save_idxs=[1,2,3])

## Numeric and analytic results
f = Figure(fontsize=18)
ax = Axis3(f[1, 1],
   title = "Curvature Drift",
   xlabel = "x [m]",
   ylabel = "y [m]",
   zlabel = "z [m]",
   aspect = :data,
   azimuth = 0.3π,
)

gc_plot(x, y, z, vx, vy, vz) = (gc(SA[x, y, z, vx, vy, vz])...,)

lines!(ax, sol, idxs=(1, 2, 3), color=Makie.wong_colors()[1])
lines!(ax, sol, idxs=(gc_plot, 1, 2, 3, 4, 5, 6), color=Makie.wong_colors()[2])
lines!(ax, sol_gc, idxs=(1, 2, 3), color=Makie.wong_colors()[3])

f = DisplayAs.PNG(f) #hide