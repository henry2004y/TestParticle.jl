# # Hybrid GC-FO Solver
#
# This example demonstrates the adaptive hybrid solver that dynamically switches
# between full orbit (FO) and guiding center (GC) tracing based on the local
# adiabaticity parameter `ε = ρ_L / R_c`.
#
# We use a magnetic bottle field where the adiabaticity varies spatially: away
# from the mirror points the field lines are strongly curved (`ε` large,
# non-adiabatic → FO), while closer to the mirror points the curvature drops
# (`ε` small, adiabatic → GC).

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
using BorisPushers
using TestParticle: get_adiabaticity
import TestParticle as TP
using LinearAlgebra, Random, Printf, Markdown
using Chairmarks, Statistics
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

const B0 = 1.0e-4   # [T]
const α = 1.0e-2    # [m⁻²]

function bottle_B(x, t)
    Bz = B0 * (1 + α * x[3]^2)
    Bx = -B0 * α * x[1] * x[3]
    By = -B0 * α * x[2] * x[3]
    return SA[Bx, By, Bz]
end

m = TP.mᵢ
q = TP.qᵢ
q2m = q / m;

# ## Step 1: Full Orbit Reference Trace
#
# First we trace the proton using the standard ODE solver to confirm
# it is trapped and bounces inside the magnetic bottle.
#
# The proton is initialized at the midplane with a perpendicular velocity
# ``v_{\perp} = 5\times10^4~\mathrm{m/s}`` in the ``x``-direction and a parallel
# velocity ``v_{\parallel} = 1\times10^5~\mathrm{m/s}`` in the ``z``-direction.
# This gives a significant perpendicular component needed for the mirror force
# to reflect the particle.
#
# The starting point is put off the magnetic axis. On the axis the field is
# exactly uniform, so its curvature — and hence `ε` — vanishes there;
# launching at ``x = 10~\mathrm{m}`` (about two gyroradii) puts the particle in
# the region where the field lines actually bend.

x0 = SA[10.0, 0.0, 0.0]
v_perp = 5.0e4  # [m/s]
v_par = 1.0e5   # [m/s]
v0 = SA[v_perp, 0.0, v_par]
u0 = vcat(x0, v0)

Ω = abs(q2m) * B0
T_gyro = 2π / Ω

## Trace long enough to see several bounces
tspan = (0.0, 30 * T_gyro)

param_fo = prepare(; B = bottle_B, species = Proton)
prob_fo = ODEProblem(trace, u0, tspan, param_fo)
sol_fo = solve(prob_fo, Vern6())

## Verify trapping: z should oscillate, not diverge
z_fo = [u[3] for u in sol_fo.u]
@assert maximum(abs, z_fo) < 1.0e4 "Particle escaped the bottle!"

# ## Step 2: Guiding Center Reference Trace
#
# For comparison, we also trace with the guiding center equations.
# Note that [`trace_gc!`](@ref) traces the guiding-center position ``\mathbf{X}``
# (a smooth curve) under the assumption that the magnetic moment ``\mu`` is strictly
# conserved, whereas full-orbit tracing resolves the fast gyration of the particle.

bottle_B_static(x) = bottle_B(x, 0.0)

stateinit_gc, param_gc = TP.prepare_gc(
    u0; B = bottle_B_static, species = Proton
)
prob_gc = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
sol_gc = solve(prob_gc, Vern6());

# ## Step 3: Hybrid Solver
#
# The adaptive hybrid solver switches between FO and GC based on the
# adiabaticity parameter and a threshold.

## The classical adiabaticity criterion: ε = ρ_L / R_c < 0.1
threshold = 0.1

alg = AdaptiveHybrid(;
    threshold,
    dtmax = T_gyro,
    dtmin = 1.0e-4 * T_gyro,
    maxiters = 500_000,
    check_interval = 10,
)

# Both modes adapt their step to the local field: the full-orbit mode follows
# the gyroperiod, the guiding-center mode uses an embedded RK45 estimate. The
# adiabaticity is checked after every `check_interval` accepted steps *or*
# every `check_dt` of simulated time, whichever comes first — a guiding-center
# step can otherwise span many gyroperiods and enter a non-adiabatic region
# unnoticed. `check_dt` defaults to `100 * dtmax`; lower it where `ε` varies
# faster than that.

## Set verbose = true to see the dynamic switching
p = param_fo
prob_hybrid = TraceProblem(u0, tspan, param_fo)
sol = solve(prob_hybrid, alg; verbose = false, seed = 1234).u[1];

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
# All three trajectories are sampled on the same uniform time grid before plotting.
# `sol_fo` and `sol_gc` are Vern6 ODE solutions with high-order dense output; the
# native hybrid solution interpolates linearly between its saved nodes, so its
# curve keeps the straight segments of its much coarser output. Sampling more
# densely does not add information: that coarseness is the point of the hybrid
# solver — it saves far fewer points — not a measure of accuracy. Shared axis
# limits are taken from the drawn curves, so no trajectory is clipped.

const npts_plot = 5000

"Sample `sol` and its adiabaticity on a uniform time grid for plotting."
function trajectory_points(sol, ε_vals; npts = npts_plot)
    t_new = range(tspan[1], tspan[2]; length = npts)
    xs = Vector{Float64}(undef, npts)
    ys = Vector{Float64}(undef, npts)
    zs = Vector{Float64}(undef, npts)
    ε_plot = Vector{Float64}(undef, npts)
    for (i, t) in enumerate(t_new)
        u = sol(t)
        xs[i], ys[i], zs[i] = u[1], u[2], u[3]
        idx = clamp(searchsortedlast(sol.t, t), 1, length(sol.t) - 1)
        frac = (t - sol.t[idx]) / (sol.t[idx + 1] - sol.t[idx])
        ε_plot[i] = ε_vals[idx] + frac * (ε_vals[idx + 1] - ε_vals[idx])
    end
    return (; xs, ys, zs, ε = ε_plot)
end

"Axis limits shared by a set of sampled trajectories."
function shared_limits(traces)
    lo = fill(Inf, 3)
    hi = fill(-Inf, 3)
    for tr in traces, (k, c) in enumerate((tr.xs, tr.ys, tr.zs))
        lo[k] = min(lo[k], minimum(c))
        hi[k] = max(hi[k], maximum(c))
    end
    return ((lo[1], hi[1]), (lo[2], hi[2]), (lo[3], hi[3]))
end

## Plot a sampled trajectory colored by ε onto an existing Axis3
function plot_trajectory!(ax, tr)
    ε_log = log10.(clamp.(tr.ε, ε_clamp_lo, Inf))
    lines!(
        ax, tr.xs, tr.ys, tr.zs;
        color = ε_log, colormap = :turbo, colorrange = clims, linewidth = 1.0
    )
    scatter!(
        ax, [Point3f(tr.xs[1], tr.ys[1], tr.zs[1])];
        color = :purple, markersize = 12
    )
    return scatter!(
        ax, [Point3f(tr.xs[end], tr.ys[end], tr.zs[end])];
        color = :red, marker = :rect, markersize = 12
    )
end

f = Figure(; size = (1400, 900), fontsize = 20)

traces = (
    trajectory_points(sol_fo, ε_fo),
    trajectory_points(sol_gc, ε_gc),
    trajectory_points(sol, ε_hybrid),
)
lims = shared_limits(traces)

## Top row: three interpolated 3D trajectories + shared colorbar
ax_fo = Axis3(
    f[1, 1],
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    title = "Full Orbit (Particle)", aspect = :data,
    limits = lims,
)
plot_trajectory!(ax_fo, traces[1])

ax_gc = Axis3(
    f[1, 2],
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    title = "Guiding Center (GC Position)", aspect = :data,
    limits = lims,
)
plot_trajectory!(ax_gc, traces[2])

ax_hyb = Axis3(
    f[1, 3],
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    title = "Hybrid Mode (Particle)", aspect = :data,
    limits = lims,
)
plot_trajectory!(ax_hyb, traces[3])

Colorbar(f[1, 4]; colormap = :turbo, limits = clims, label = L"\log_{10}(\epsilon)")

## Bottom row: time series of adiabaticity with solver mode
t_norm = sol.t ./ T_gyro

ax_ts = Axis(
    f[2, 1:4],
    xlabel = L"t / T_\text{gyro}",
    ylabel = L"\epsilon = \rho_L / R_c",
    yscale = log10,
    title = "Adiabaticity & Solver Mode",
    limits = (
        (minimum(t_norm), maximum(t_norm)),
        (ε_clamp_lo, 10^ceil(log10(maximum(ε_hybrid)))),
    ),
)

## Shade FO and GC regions using the solver's recorded modes
let
    t_adia = sol.stats.adiabaticity.t ./ T_gyro
    mode_adia = sol.stats.adiabaticity.mode
    i_region = 1
    while i_region <= length(t_adia)
        mode_fo = mode_adia[i_region] === :FO
        j_region = i_region
        while j_region < length(t_adia) && (mode_adia[j_region + 1] === :FO) == mode_fo
            j_region += 1
        end
        t_lo = t_adia[i_region]
        t_hi = j_region < length(t_adia) ? t_adia[j_region + 1] : t_norm[end]
        c = mode_fo ? (:red, 0.2) : (:blue, 0.2)
        vspan!(ax_ts, t_lo, t_hi; color = c)
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

# ## Performance Comparison
#
# Finally, we compare the execution time, memory and output size of the three
# solvers. Pointwise position differences are a poor accuracy measure in this
# bottle: the gyrophase and the bounce phase are exponentially sensitive, and the
# guiding-center equations describe a different curve from the particle. The drift
# of the kinetic energy is the robust numerical accuracy measure here.

"Largest relative energy drift over the trajectory."
function energy_drift(sol)
    e = map(sol.u) do u
        if length(u) == 6
            return 0.5 * m * sum(abs2, u[SA[4, 5, 6]])
        end
        x = SA[u[1], u[2], u[3]]
        return 0.5 * m * u[4]^2 + param_gc.μ * norm(bottle_B_static(x))
    end
    return (maximum(e) - minimum(e)) / e[1]
end

b_fo = @be solve(prob_fo, Vern6())
b_gc = @be solve(prob_gc, Vern6())
b_hy = @be solve(prob_hybrid, alg; verbose = false)

io = IOBuffer() #hide
println(io, "| Solver | Time | Allocations | Saved points | Energy drift |") #hide
println(io, "| :--- | :--- | :--- | ---: | ---: |") #hide
Printf.@printf(
    io, "| Full Orbit | %.2f μs | %.2f KiB | %d | %.1e |\n",
    median(b_fo).time * 1.0e6, median(b_fo).bytes / 1024,
    length(sol_fo.t), energy_drift(sol_fo)
) #hide
Printf.@printf(
    io, "| Guiding Center | %.2f μs | %.2f KiB | %d | %.1e |\n",
    median(b_gc).time * 1.0e6, median(b_gc).bytes / 1024,
    length(sol_gc.t), energy_drift(sol_gc)
) #hide
Printf.@printf(
    io, "| Hybrid | %.2f μs | %.2f KiB | %d | %.1e |\n",
    median(b_hy).time * 1.0e6, median(b_hy).bytes / 1024,
    length(sol.t), energy_drift(sol)
) #hide
Markdown.parse(String(take!(io))) #hide

# ## AdaptiveHybrid: Adiabaticity Check Modes
#
# [`AdaptiveHybrid`](@ref) lets the user pick which adiabaticity criterion drives
# the GC ↔ FO decision, via the `adiabaticity` keyword:
# - `:curvature` (default) → `ε_curv = ρ_L / R_c` (curvature drift; reproduces
#   the legacy behaviour),
# - `:gradB`     → `ε_gradB = ρ_L / L_B`, with `L_B = |B| / |∇B|`,
# - `:both`      → OR of the two criteria: switch to full orbit whenever *either*
#   `ε_curv ≥ α` *or* `ε_gradB ≥ α` (equivalently `max(ε_curv, ε_gradB) ≥ α`),
# - `:jacobian`  → `ε_jac = ρ_L · ‖JB‖_F / |B|`, all-in-one criterion, which also
#   captures torsion and shear.
#
# The criteria solve the same problem; they differ only in *when* the solver drops
# into the full-orbit mode, so the trajectories stay close while the time spent
# in each mode changes.

mode_threshold = 0.1
mode_common = (;
    threshold = mode_threshold, dtmax = T_gyro,
    dtmin = 1.0e-4 * T_gyro, maxiters = 500_000, check_interval = 10,
)
alg_curv = AdaptiveHybrid(; mode_common..., adiabaticity = :curvature)
alg_gradB = AdaptiveHybrid(; mode_common..., adiabaticity = :gradB)
alg_both = AdaptiveHybrid(; mode_common..., adiabaticity = :both)

prob_mode = TraceProblem(u0, tspan, p)
sol_curv = solve(prob_mode, alg_curv; verbose = false, seed = 1234).u[1]
sol_gradB = solve(prob_mode, alg_gradB; verbose = false, seed = 1234).u[1]
sol_both = solve(prob_mode, alg_both; verbose = false, seed = 1234).u[1];

# RMS position difference vs the full-orbit reference (dominated by gyrophase and
# bounce phase here), and the fraction of simulated time spent in the FO mode.
function _adia_rms(sol)
    errs = Float64[]
    for (t, u) in zip(sol.t, sol.u)
        u_ref = sol_fo(t)
        push!(errs, norm(u[1:3] - u_ref[1:3]))
    end
    return sqrt(mean(errs .^ 2))
end
function _adia_fo_frac(sol)
    ad = sol.stats.adiabaticity
    isempty(ad.t) && return 0.0
    t_end = sol.t[end]
    t_fo = 0.0
    for i in eachindex(ad.t)
        t_hi = i < length(ad.t) ? ad.t[i + 1] : t_end
        ad.mode[i] === :FO && (t_fo += t_hi - ad.t[i])
    end
    return t_fo / (t_end - ad.t[1])
end

ε_col_curv = get_adiabaticity(sol_curv)
ε_col_gradB = get_adiabaticity(sol_gradB)
ε_col_both = get_adiabaticity(sol_both)

io2 = IOBuffer() #hide
println(io2, "### Result vs full-orbit reference") #hide
println(io2, "") #hide
println(io2, "| Mode | FO time fraction | Saved points | RMS pos. diff [m] |") #hide
println(io2, "| :--- | :--- | :--- | :--- |") #hide
Printf.@printf(
    io2, "| `:curvature` | %.2f | %d | %.2e |\n",
    _adia_fo_frac(sol_curv), length(sol_curv.t), _adia_rms(sol_curv)
) #hide
Printf.@printf(
    io2, "| `:gradB` | %.2f | %d | %.2e |\n",
    _adia_fo_frac(sol_gradB), length(sol_gradB.t), _adia_rms(sol_gradB)
) #hide
Printf.@printf(
    io2, "| `:both` | %.2f | %d | %.2e |\n",
    _adia_fo_frac(sol_both), length(sol_both.t), _adia_rms(sol_both)
) #hide
Markdown.parse(String(take!(io2))) #hide

# ### Trajectories under the three modes
#
# The three modes trace nearly the same path (all stay close to the full-orbit
# reference); they differ in *when* the solver drops into the full-orbit mode.
# As above, every curve is interpolated onto the same time grid and the shared
# limits come from those curves, so the three panels are drawn on one scale.

f_modes = Figure(; size = (1400, 460), fontsize = 20)

traces_modes = (
    trajectory_points(sol_curv, ε_col_curv),
    trajectory_points(sol_gradB, ε_col_gradB),
    trajectory_points(sol_both, ε_col_both),
)
lms = shared_limits(traces_modes)

ax_c = Axis3(
    f_modes[1, 1], title = ":curvature", aspect = :data, limits = lms,
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]"
)
plot_trajectory!(ax_c, traces_modes[1])

ax_g = Axis3(
    f_modes[1, 2], title = ":gradB", aspect = :data, limits = lms,
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]"
)
plot_trajectory!(ax_g, traces_modes[2])

ax_b = Axis3(
    f_modes[1, 3], title = ":both", aspect = :data, limits = lms,
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]"
)
plot_trajectory!(ax_b, traces_modes[3])

f_modes = DisplayAs.PNG(f_modes) #hide

# ### Why the modes differ: the selected adiabaticity
#
# The diagnostics store *only* the adiabaticity value used by the selected
# mode (one scalar per check point), not the full component vector. For
# `:curvature` and `:gradB` that is `ε_curv` and `ε_gradB`; for `:both` it is
# `max(ε_curv, ε_gradB)`. Plotting the stored value for each mode shows that
# `:both` switches at least as often as either single criterion alone.

function _adia_traces(sol)
    t = sol.stats.adiabaticity.t ./ T_gyro
    ε = sol.stats.adiabaticity.components
    mode = sol.stats.adiabaticity.mode
    return t, ε, mode
end
t_curv, ε_curv, mode_curv = _adia_traces(sol_curv)
t_gradB, ε_gradB, mode_gradB = _adia_traces(sol_gradB)
t_both, ε_both, mode_both = _adia_traces(sol_both)

f_comp = Figure(; size = (1000, 480), fontsize = 20)
ax_comp = Axis(
    f_comp[1, 1],
    xlabel = L"t / T_\text{gyro}",
    ylabel = L"\epsilon",
    yscale = log10,
    title = "Selected adiabaticity by mode",
    limits = (
        (0.0, 30.0),
        (ε_clamp_lo, 10.0),
    ),
)

## Shade full-orbit (red) vs guiding-center (blue) regions for each mode.
function _shade!(ax, t, mode)
    i = 1
    while i <= length(t)
        m_fo = mode[i] === :FO
        j = i
        while j < length(t) && mode[j + 1] === mode[i]
            j += 1
        end
        c = m_fo ? (:red, 0.18) : (:blue, 0.18)
        vspan!(ax, t[i], t[j]; color = c)
        i = j + 1
    end
    return
end
_shade!(ax_comp, t_curv, mode_curv)
_shade!(ax_comp, t_gradB, mode_gradB)
_shade!(ax_comp, t_both, mode_both)

lines!(ax_comp, t_curv, ε_curv; color = :green, label = ":curvature (ε_curv)")
lines!(ax_comp, t_gradB, ε_gradB; color = :orange, label = ":gradB (ε_gradB)")
lines!(ax_comp, t_both, ε_both; color = :black, linewidth = 1.6, label = ":both (max)")
hlines!(
    ax_comp, [mode_threshold];
    color = :gray50, linestyle = :dash, linewidth = 1.5,
    label = "threshold = $(round(mode_threshold; sigdigits = 2))"
)
axislegend(ax_comp; position = :rt)
f_comp = DisplayAs.PNG(f_comp) #hide

# ### Performance
#
# All three modes share the same solver core; `:curvature` carries only a few KiB
# of extra diagnostics. The `:gradB` / `:both` modes can even run *faster* in this
# setup because they spend more time in the cheap Boris full-orbit integrator, at
# the cost of higher memory (more saved points).

b_curv = @be solve($prob_mode, $alg_curv; verbose = false, seed = 1234)
b_gradB = @be solve($prob_mode, $alg_gradB; verbose = false, seed = 1234)
b_both = @be solve($prob_mode, $alg_both; verbose = false, seed = 1234)

io3 = IOBuffer() #hide
println(io3, "| Mode | Time | Allocations |") #hide
println(io3, "| :--- | :--- | :--- |") #hide
Printf.@printf(
    io3, "| `:curvature` | %.2f ms | %.2f KiB |\n",
    median(b_curv).time * 1.0e3, median(b_curv).bytes / 1024
) #hide
Printf.@printf(
    io3, "| `:gradB` | %.2f ms | %.2f KiB |\n",
    median(b_gradB).time * 1.0e3, median(b_gradB).bytes / 1024
) #hide
Printf.@printf(
    io3, "| `:both` | %.2f ms | %.2f KiB |\n",
    median(b_both).time * 1.0e3, median(b_both).bytes / 1024
) #hide
Markdown.parse(String(take!(io3))) #hide
