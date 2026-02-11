# # Hybrid GC-FO Solver
#
# This example demonstrates the adaptive hybrid solver that dynamically switches
# between full orbit (FO) and guiding center (GC) tracing based on the local
# adiabaticity parameter `ε = ρ_L / R_c`.
#
# We use a magnetic bottle field where the adiabaticity varies spatially:
# near the midplane `ε` is large (non-adiabatic → FO), and near the mirror
# points `ε` is small (adiabatic → GC).

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
import TestParticle as TP
using LinearAlgebra, Random
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Magnetic Bottle Configuration
#
# The magnetic field satisfies `∇·B = 0`:
# ```math
# B_z = B_0 (1 + α z^2), \quad
# B_x = -B_0 α x z, \quad
# B_y = -B_0 α y z
# ```
# A larger `α` produces stronger curvature at the midplane, which
# raises the adiabaticity parameter `ε` there.

B0 = 1.0e-4   # [T]
α = 1.0e-2    # [m⁻²]

function bottle_B(x, t)
    Bz = B0 * (1 + α * x[3]^2)
    Bx = -B0 * α * x[1] * x[3]
    By = -B0 * α * x[2] * x[3]
    return SA[Bx, By, Bz]
end

B_field = TP.Field(bottle_B)
E_field = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

m = TP.mᵢ
q = TP.qᵢ
q2m = q / m;

# ## Step 1: Full Orbit Reference Trace
#
# First we trace the proton using the standard ODE solver to confirm
# it is trapped and bounces inside the magnetic bottle.

x0 = SA[0.0, 0.0, 0.0]
v_perp = 5.0e4  # [m/s]
v_par = 1.0e5   # [m/s]
v0 = SA[v_perp, 0.0, v_par]
u0 = vcat(x0, v0)

Ω = abs(q2m) * B0
T_gyro = 2π / Ω

## Trace long enough to see several bounces
tspan = (0.0, 30 * T_gyro)

param_fo = prepare(E_field, B_field; species = Proton)
prob_fo = ODEProblem(trace, u0, tspan, param_fo)
sol_fo = solve(prob_fo, Vern9())

## Verify trapping: z should oscillate, not diverge
z_fo = [u[3] for u in sol_fo.u]
@assert maximum(abs, z_fo) < 1.0e4 "Particle escaped the bottle!"

# ## Step 2: Guiding Center Reference Trace
#
# For comparison, we also trace with the guiding center equations.

bottle_B_static(x) = bottle_B(x, 0.0)
bottle_E_static(x) = SA[0.0, 0.0, 0.0]

stateinit_gc, param_gc = TP.prepare_gc(
    u0, bottle_E_static, bottle_B_static; species = Proton
)
prob_gc = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
sol_gc = solve(prob_gc, Vern9());

# ## Step 3: Hybrid Solver
#
# The adaptive hybrid solver switches between FO and GC based on the
# adiabaticity parameter and a threshold.

## The classical adiabaticity criterion: ε = ρ_L / R_c < 0.1
threshold = 0.1

p = (q2m, m, E_field, B_field, TP.ZeroField())

alg = AdaptiveHybrid(;
    threshold,
    dtmax = T_gyro,
    dtmin = 1.0e-4 * T_gyro,
    maxiters = 500_000,
)

Random.seed!(1234)
## Set verbose = true to see the dynamic switching
sol = TP.solve(TraceHybridProblem(u0, tspan, p), alg; verbose = false)[1];

# ## Step 4: Compute Adiabaticity
#
# The adiabaticity parameter `ε = ρ_L / R_c` measures how well the
# guiding center approximation holds. When `ε < threshold`, GC is
# valid; when `ε ≥ threshold`, we need full orbit tracing.

ε_fo = get_adiabaticity(sol_fo)
ε_gc = get_adiabaticity(sol_gc)
ε_hybrid = get_adiabaticity(sol)

## Common colorbar range across all solvers
ε_clamp_lo = 1.0e-3
clims = extrema(log10.(clamp.(vcat(ε_fo, ε_gc, ε_hybrid), ε_clamp_lo, Inf)))

# ## Visualization
#
# We define a helper to plot a 3D trajectory colored by ε, with
# start/end markers and a colorbar.

## Plot a 3D trajectory colored by ε onto an existing Axis3
## When `npts` is set, the solution is interpolated for a smoother curve.
function plot_trajectory!(ax, sol, ε_vals; npts = nothing)
    if isnothing(npts)
        xs = [u[1] for u in sol.u]
        ys = [u[2] for u in sol.u]
        zs = [u[3] for u in sol.u]
        ε_plot = ε_vals
    else
        t_new = range(sol.t[1], sol.t[end]; length = npts)
        xs = Vector{Float64}(undef, npts)
        ys = Vector{Float64}(undef, npts)
        zs = Vector{Float64}(undef, npts)
        ε_plot = Vector{Float64}(undef, npts)
        for (i, t) in enumerate(t_new)
            u = sol(t)
            xs[i], ys[i], zs[i] = u[1], u[2], u[3]
            ## Linearly interpolate ε to the new time grid
            idx = clamp(searchsortedlast(sol.t, t), 1, length(sol.t) - 1)
            frac = (t - sol.t[idx]) / (sol.t[idx + 1] - sol.t[idx])
            ε_plot[i] = ε_vals[idx] + frac * (ε_vals[idx + 1] - ε_vals[idx])
        end
    end
    ε_log = log10.(clamp.(ε_plot, ε_clamp_lo, Inf))
    lines!(
        ax, xs, ys, zs;
        color = ε_log, colormap = :turbo, colorrange = clims, linewidth = 1.0
    )
    scatter!(ax, [Point3f(xs[1], ys[1], zs[1])]; color = :purple, markersize = 12)
    return scatter!(
        ax, [Point3f(xs[end], ys[end], zs[end])];
        color = :red, marker = :rect, markersize = 12
    )
end

## Create a standalone figure with a 3D trajectory colored by ε
function plot_trajectory(sol, ε_vals, title_str; figsize = (700, 500), npts = nothing)
    f = Figure(; size = figsize, fontsize = 18)
    ax = Axis3(
        f[1, 1],
        xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
        title = title_str, aspect = :data,
    )
    plot_trajectory!(ax, sol, ε_vals; npts)
    Colorbar(f[1, 2]; colormap = :turbo, limits = clims, label = L"\log_{10}(\epsilon)")
    return f
end

# ## Visualization

f = Figure(; size = (1400, 900), fontsize = 18)

## Top row: three 3D trajectories + shared colorbar
ax_fo = Axis3(
    f[1, 1],
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    title = "Full Orbit", aspect = :data,
)
plot_trajectory!(ax_fo, sol_fo, ε_fo; npts = 5000)

ax_gc = Axis3(
    f[1, 2],
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    title = "Guiding Center", aspect = :data,
)
plot_trajectory!(ax_gc, sol_gc, ε_gc; npts = 5000)

ax_hyb = Axis3(
    f[1, 3],
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    title = "Hybrid Mode", aspect = :data,
)
plot_trajectory!(ax_hyb, sol, ε_hybrid)

Colorbar(f[1, 4]; colormap = :turbo, limits = clims, label = L"\log_{10}(\epsilon)")

## Bottom row: time series of adiabaticity with solver mode
ax_ts = Axis(
    f[2, 1:4],
    xlabel = L"t / T_\text{gyro}",
    ylabel = L"\epsilon = \rho_L / R_c",
    yscale = log10,
    title = "Adiabaticity & Solver Mode",
    limits = (nothing, (ε_clamp_lo, 10^ceil(log10(maximum(ε_hybrid))))),
)

t_norm = sol.t ./ T_gyro

## Shade FO and GC regions
is_fo = ε_hybrid .>= threshold
let i_region = 1
    while i_region <= length(t_norm)
        mode_fo = is_fo[i_region]
        j_region = i_region
        while j_region < length(t_norm) && is_fo[j_region + 1] == mode_fo
            j_region += 1
        end
        t_lo = t_norm[i_region]
        t_hi = t_norm[j_region]
        if mode_fo
            vspan!(ax_ts, t_lo, t_hi; color = (:red, 0.2))
        else
            vspan!(ax_ts, t_lo, t_hi; color = (:blue, 0.2))
        end
        i_region = j_region + 1
    end
end

lines!(ax_ts, t_norm, ε_hybrid; color = :black, linewidth = 1.0)
hlines!(
    ax_ts, [threshold];
    color = :gray50, linestyle = :dash, linewidth = 1.5,
    label = "threshold = $(round(threshold; sigdigits = 2))"
)

## Legend entries for solver modes
poly!(
    ax_ts, Point2f[(NaN, NaN)];
    color = (:red, 0.4), strokewidth = 0, label = "Full Orbit"
)
poly!(
    ax_ts, Point2f[(NaN, NaN)];
    color = (:blue, 0.4), strokewidth = 0,
    label = "Guiding Center"
)
Legend(f[3, 1:4], ax_ts; orientation = :horizontal)

rowsize!(f.layout, 1, Relative(0.55))

f = DisplayAs.PNG(f) #hide
