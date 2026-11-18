# # Hybrid GC-FO Solver
#
# This example demonstrates [`AdaptiveHybrid`](@ref), which dynamically switches between
# full-orbit (FO) and guiding-center (GC) tracing based on a local adiabaticity parameter
# such as `ε = ρ_L / R_c`.
#
# We use a divergence-free curved magnetic loop where the field strength varies by nearly
# two orders of magnitude along circular field lines: in the strong-field mirror legs `ε`
# is small (adiabatic → GC), whereas near the weak-field loop apex `ε` exceeds the
# switching threshold (non-adiabatic → FO).

import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using BorisPushers
using TestParticle: get_adiabaticity
import TestParticle as TP
using LinearAlgebra, Printf, Markdown
using Chairmarks, Statistics
using CairoMakie
CairoMakie.activate!(; type = "png") #hide

# ## Curved Magnetic Loop Configuration
#
# Derived from the 2-D line-dipole flux function
# ```math
# \psi(x, z) = -\frac{M X}{X^2 + z^2}, \qquad X = x + x_0,
# ```
# the magnetic field ``\mathbf{B} = \nabla \times (\psi \hat{\mathbf{y}})`` satisfies
# ``\nabla \cdot \mathbf{B} = 0`` identically:
# ```math
# B_x = -\frac{2 M X z}{r^4}, \quad
# B_y = 0, \quad
# B_z = \frac{M (X^2 - z^2)}{r^4}, \qquad r^2 = X^2 + z^2.
# ```
# Every field line ``X^2 + z^2 = D X`` is a circle of constant curvature radius
# ``R_c = D / 2`` whose field magnitude ``|\mathbf{B}| = M / r^2`` grows strongly toward
# the footpoints near ``(-x_0, 0)`` and weakens at the apex ``(D - x_0, 0)``.

const M_LOOP = 1.6e-5   # [T m²]; gives |B| = 1 μT at the D = 4 m apex
const X_OFFSET = 3.0    # [m]; line-dipole origin at (x = -3 m, z = 0)
const D_LINE = 4.0      # [m]; field-line circle diameter

function curved_B(x, t = 0.0)
    X = x[1] + X_OFFSET
    z = x[3]
    r2 = X * X + z * z
    r4 = r2 * r2
    return SA[
        -2 * M_LOOP * X * z / r4,
        0.0,
        M_LOOP * (X * X - z * z) / r4,
    ]
end

m = TP.mᵢ
q = TP.qᵢ
q2m = q / m;

# ## Step 1: Full-Orbit Reference Trace
#
# We launch a proton in the southern strong-field leg (`r_start = 0.75 m` on the
# `D = 4 m` field line, where `|B|` is `28×` stronger than at the apex) with speed
# `42 m/s` and pitch angle `65°`. The proton bounces between conjugate mirror points while
# curvature and `∇B` drifts carry it along `-y`.

r_start = 0.75
v_tot = 42.0
pitch_deg = 65.0

X0 = r_start^2 / D_LINE
z0 = -sqrt(r_start^2 - X0^2)
x0 = SA[X0 - X_OFFSET, 0.0, z0]
B0_vec = curved_B(x0)
B0_mag = norm(B0_vec)
b0 = B0_vec / B0_mag
α0 = deg2rad(pitch_deg)
v0 = v_tot * (cos(α0) * b0 + sin(α0) * SA[0.0, 1.0, 0.0])
u0 = vcat(x0, v0)

B_apex = M_LOOP / D_LINE^2
T_gyro_0 = 2π / (abs(q2m) * B0_mag)
T_gyro_apex = 2π / (abs(q2m) * B_apex)
tspan = (0.0, 1.2)

param_fo = prepare(; B = curved_B, species = Proton)
prob_fo = TraceProblem(u0, tspan, param_fo)
sol_fo = solve(prob_fo, Vern7(); abstol = 1.0e-7, reltol = 1.0e-7, dense = true)

alg_boris = AdaptiveBoris(; safety = 0.04)
sol_boris = solve(
    prob_fo, alg_boris; dtmax = T_gyro_apex, dtmin = 1.0e-3 * T_gyro_0
);

# ## Step 2: Guiding-Center Reference Trace
#
# Standalone guiding-center integration evolves ``(\mathbf{X}, v_\parallel)`` assuming the
# first adiabatic invariant ``\mu`` remains strictly constant across the weak-field apex.

u0_gc, param_gc = prepare_gc(u0; B = (x) -> curved_B(x), species = Proton)
prob_gc = TraceGCProblem(SVector{4}(u0_gc), tspan, param_gc)
sol_gc = solve(prob_gc, Vern7(); abstol = 1.0e-7, reltol = 1.0e-7, dense = true);

# ## Step 3: Adaptive Hybrid Solver
#
# [`AdaptiveHybrid`](@ref) uses hysteresis thresholds (`threshold_gc_to_fo ≥
# threshold_fo_to_gc`) to switch cleanly from GC to FO entering the curved weak-field apex
# and back to GC in the strong-field mirror legs without chattering.

thr_gc_to_fo = 0.1
thr_fo_to_gc = 0.08

alg = AdaptiveHybrid(;
    threshold_gc_to_fo = thr_gc_to_fo,
    threshold_fo_to_gc = thr_fo_to_gc,
    dtmax = T_gyro_apex,
    dtmin = 1.0e-2 * T_gyro_0,
    safety_fo = 0.04,
    check_dt = 0.25 * T_gyro_apex,
    check_interval = 4,
    abstol = 1.0e-7,
    reltol = 1.0e-7,
    maxiters = 500_000,
)

sol = solve(prob_fo, alg; verbose = false, seed = 42).u[1];

# ## Step 4: Visualization
#
# Comparing the 3-D trajectories alongside the poloidal `(x–z)` projection and the
# recorded mode timeline shows that pure GC misses the non-adiabatic mirror shift across
# the apex, whereas `AdaptiveHybrid` tracks the full-orbit trajectory while stepping
# rapidly through the strong-field mirror bounces.

ε_fo = get_adiabaticity(sol_fo)
ε_gc = get_adiabaticity(sol_gc)
ε_hyb = get_adiabaticity(sol)

const ε_clamp_lo = 3.0e-3
clims = (log10(ε_clamp_lo), log10(0.35))

function cubic_eval(t::AbstractVector{Float64}, y::AbstractVector{Float64}, ts)
    n = length(t)
    h = diff(t)
    dl = h[1:(n - 2)]
    du = h[2:(n - 1)]
    d = 2.0 .* (dl .+ du)
    rhs = 6.0 .* ((y[3:n] .- y[2:(n - 1)]) ./ du .- (y[2:(n - 1)] .- y[1:(n - 2)]) ./ dl)
    m_c = vcat(0.0, SymTridiagonal(d, h[2:(n - 2)]) \ rhs, 0.0)
    out = Vector{Float64}(undef, length(ts))
    for (k, tq) in enumerate(ts)
        i = clamp(searchsortedlast(t, tq), 1, n - 1)
        hi = h[i]
        a = (t[i + 1] - tq) / hi
        b = (tq - t[i]) / hi
        out[k] = a * y[i] + b * y[i + 1] +
            ((a^3 - a) * m_c[i] + (b^3 - b) * m_c[i + 1]) * (hi^2 / 6.0)
    end
    return out
end

function sample_dense(s, eps_vals, npts)
    ts = range(s.t[1], s.t[end]; length = npts)
    es = cubic_eval(s.t, eps_vals, ts)
    xs = Vector{Float64}(undef, npts)
    ys = Vector{Float64}(undef, npts)
    zs = Vector{Float64}(undef, npts)
    for (i, t) in enumerate(ts)
        u = s(t)
        xs[i], ys[i], zs[i] = u[1], u[2], u[3]
    end
    return xs, ys, zs, es
end

function sample_spline(s, eps_vals, npts)
    keep = [i == 1 || s.t[i] > s.t[i - 1] for i in eachindex(s.t)]
    t_u = s.t[keep]
    u_u = s.u[keep]
    e_u = eps_vals[keep]
    ts = range(t_u[1], t_u[end]; length = npts)
    return (
        cubic_eval(t_u, [u[1] for u in u_u], ts),
        cubic_eval(t_u, [u[2] for u in u_u], ts),
        cubic_eval(t_u, [u[3] for u in u_u], ts),
        cubic_eval(t_u, e_u, ts),
    )
end

function shade_modes!(ax, t, mode; alpha = 0.18)
    for i in 1:(length(t) - 1)
        t[i + 1] > t[i] || continue
        c = mode[i] === :FO ? (:red, alpha) : (:blue, alpha)
        vspan!(ax, t[i], t[i + 1]; color = c)
    end
    return nothing
end

all_u = vcat(sol_fo.u, sol_gc.u, sol.u)
pad(lu, f = 0.06) = (lu[1] - f * (lu[2] - lu[1]), lu[2] + f * (lu[2] - lu[1]))
lims3d = (
    pad(extrema(u[1] for u in all_u)),
    pad(extrema(u[2] for u in all_u)),
    pad(extrema(u[3] for u in all_u)),
)

f = Figure(; size = (1400, 880), fontsize = 20)
gl_top = f[1, 1] = GridLayout()
gl_bot = f[2, 1] = GridLayout()

xs_fo_d, ys_fo_d, zs_fo_d, es_fo_d = sample_dense(sol_fo, ε_fo, 10_000)
xs_gc_d, ys_gc_d, zs_gc_d, es_gc_d = sample_dense(sol_gc, ε_gc, 4_000)
xs_hyb, ys_hyb, zs_hyb, es_hyb = sample_spline(sol, ε_hyb, 8_000)

for (col, (ttl, xs, ys, zs, es)) in enumerate(
        (
            ("(a) Full Orbit (Vern7)", xs_fo_d, ys_fo_d, zs_fo_d, es_fo_d),
            ("(b) Guiding Center (Vern7)", xs_gc_d, ys_gc_d, zs_gc_d, es_gc_d),
            ("(c) Adaptive Hybrid", xs_hyb, ys_hyb, zs_hyb, es_hyb),
        )
    )
    ax = Axis3(
        gl_top[1, col];
        title = ttl, xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
        limits = lims3d, azimuth = -0.38π, elevation = 0.16π,
        aspect = (1.12, 1.08, 0.92),
        protrusions = (24, 34, 26, 16),
    )
    lines!(
        ax, xs, ys, zs;
        color = log10.(clamp.(es, ε_clamp_lo, 10.0)),
        colormap = :turbo, colorrange = clims, linewidth = 1.6,
    )
    scatter!(ax, [Point3f(xs[1], ys[1], zs[1])]; color = :purple, markersize = 12)
    scatter!(
        ax, [Point3f(xs[end], ys[end], zs[end])];
        color = :red, marker = :rect, markersize = 12,
    )
end

Colorbar(gl_top[1, 4]; colormap = :turbo, limits = clims, label = L"\log_{10}(\varepsilon)")
colgap!(gl_top, 3, 48)

xz_xlim = pad(extrema(vcat(xs_fo_d, xs_gc_d, xs_hyb, [-X_OFFSET, 4.9 - X_OFFSET])), 0.05)
xz_zlim = pad(extrema(vcat(zs_fo_d, zs_gc_d, zs_hyb, [-2.45, 2.45])), 0.05)
ax_xz = Axis(
    gl_bot[1, 1];
    title = "(d) Poloidal (x–z) projection",
    xlabel = "x [m]", ylabel = "z [m]",
    aspect = DataAspect(),
    limits = (xz_xlim, xz_zlim),
)
for D_fl in (2.2, 3.1, 4.0, 4.9)
    φs = range(0.14, 2π - 0.14; length = 240)
    R = D_fl / 2
    lines!(
        ax_xz, (R .- R .* cos.(φs)) .- X_OFFSET, -R .* sin.(φs);
        color = (:gray50, 0.65),
        linewidth = D_fl == D_LINE ? 1.5 : 1.0,
        linestyle = D_fl == D_LINE ? :dash : :dot,
    )
end
wc = Makie.wong_colors()
xs_gc_2d, _, zs_gc_2d, _ = sample_dense(sol_gc, ε_gc, 600)
lines!(ax_xz, xs_fo_d, zs_fo_d; color = (wc[1], 0.5), linewidth = 1.4, label = "Full orbit")
lines!(
    ax_xz, xs_gc_2d, zs_gc_2d;
    color = wc[2], linestyle = :dash, linewidth = 2.5, label = "Guiding center",
)
lines!(ax_xz, xs_hyb, zs_hyb; color = wc[3], linewidth = 2.0, label = "Adaptive hybrid")
scatter!(ax_xz, [Point2f(xs_fo_d[1], zs_fo_d[1])]; color = :purple, markersize = 11)
scatter!(
    ax_xz, [Point2f(xs_fo_d[end], zs_fo_d[end])];
    color = wc[1], marker = :rect, markersize = 11,
)
scatter!(
    ax_xz, [Point2f(xs_gc_d[end], zs_gc_d[end])];
    color = wc[2], marker = :rect, markersize = 11,
)
scatter!(
    ax_xz, [Point2f(xs_hyb[end], zs_hyb[end])];
    color = wc[3], marker = :dtriangle, markersize = 12,
)
axislegend(ax_xz; position = :lb, backgroundcolor = (:white, 0.9), labelsize = 16)

t_norm = sol.stats.adiabaticity.t ./ T_gyro_apex
eps_t = sol.stats.adiabaticity.components
mode_t = sol.stats.adiabaticity.mode

ax_ts = Axis(
    gl_bot[1, 2];
    title = "(e) Local adiabaticity and active solver mode",
    xlabel = L"t / T_{\mathrm{gyro,apex}}",
    ylabel = L"\varepsilon = \rho_L / R_c",
    yscale = log10,
    limits = ((t_norm[1], t_norm[end]), (ε_clamp_lo, 1.0)),
)
shade_modes!(ax_ts, t_norm, mode_t)
lines!(ax_ts, t_norm, clamp.(eps_t, ε_clamp_lo, 1.0); color = :black, linewidth = 1.7)
hlines!(
    ax_ts, [thr_gc_to_fo];
    color = :gray35, linestyle = :dash, linewidth = 1.5,
    label = L"\alpha = 0.10\ (\mathrm{GC}\to\mathrm{FO})",
)
hlines!(
    ax_ts, [thr_fo_to_gc];
    color = :gray45, linestyle = :dot, linewidth = 1.5,
    label = L"\beta = 0.08\ (\mathrm{FO}\to\mathrm{GC})",
)
poly!(ax_ts, Point2f[(NaN, NaN)]; color = (:red, 0.4), label = "Full orbit")
poly!(ax_ts, Point2f[(NaN, NaN)]; color = (:blue, 0.4), label = "Guiding center")
axislegend(
    ax_ts; position = :rt, nbanks = 2, backgroundcolor = (:white, 0.9), labelsize = 16
)

colsize!(gl_bot, 1, Relative(0.35))
rowsize!(f.layout, 1, Relative(0.56))

f = DisplayAs.PNG(f) #hide

# ## Accuracy and Performance Comparison
#
# Because full-orbit integrators must resolve hundreds of short gyroperiods inside the
# strong-field mirror legs where the motion is deeply adiabatic, `AdaptiveHybrid` achieves
# near-GC execution speed while preserving full-orbit accuracy across the curved apex:

function energy_drift(s)
    e = map(s.u) do u
        length(u) == 6 && return 0.5 * m * sum(abs2, u[SA[4, 5, 6]])
        x = SA[u[1], u[2], u[3]]
        return 0.5 * m * u[4]^2 + param_gc.μ * norm(curved_B(x))
    end
    return (maximum(e) - minimum(e)) / e[1]
end

end_pos_err(s) = norm(s.u[end][1:3] - sol_fo.u[end][1:3])

b_fo = @be solve($prob_fo, Vern7(); abstol = 1.0e-7, reltol = 1.0e-7, dense = false)
b_bo = @be solve(
    $prob_fo, $alg_boris; dtmax = $T_gyro_apex, dtmin = $(1.0e-3 * T_gyro_0)
)
b_gc = @be solve($prob_gc, Vern7(); abstol = 1.0e-7, reltol = 1.0e-7, dense = false)
b_hy = @be solve($prob_fo, $alg; verbose = false, seed = 42)

io = IOBuffer() #hide
println( #hide
    io, #hide
    "| Solver | Time | Allocations | Saved pts | End-pos diff [m] | Energy drift |", #hide
) #hide
println(io, "| :--- | ---: | ---: | ---: | ---: | ---: |") #hide
for (label, b, s) in ( #hide
        ("Full Orbit (Vern7)", b_fo, sol_fo), #hide
        ("Full Orbit (AdaptiveBoris)", b_bo, sol_boris), #hide
        ("Guiding Center (Vern7)", b_gc, sol_gc), #hide
        ("AdaptiveHybrid", b_hy, sol), #hide
    ) #hide
    Printf.@printf( #hide
        io, "| %s | %.1f μs | %.1f KiB | %d | %.2f | %.1e |\n", #hide
        label, median(b).time * 1.0e6, median(b).bytes / 1024, #hide
        length(s.t), end_pos_err(s), energy_drift(s), #hide
    ) #hide
end #hide
Markdown.parse(String(take!(io))) #hide

# ## Adiabaticity Check Modes
#
# [`AdaptiveHybrid`](@ref) lets the user select which local criterion drives the GC ↔ FO
# decision via the `adiabaticity` keyword:
# - `:curvature` (default): `ε_curv = ρ_L / R_c`, scaling as `(r / D)²` along a circular
#   loop of radius `D / 2`;
# - `:gradB`: `ε_gradB = ρ_L / L_B` with `L_B = |B| / |∇B| = r / 2`, scaling linearly as
#   `r / D` and therefore switching earlier in the converging legs;
# - `:both`: `max(ε_curv, ε_gradB)`;
# - `:jacobian`: `ε_jac = ρ_L · ‖∇B‖_F / |B|`, capturing curvature, gradient,
#   torsion, and shear simultaneously.

mode_kw = (;
    threshold_gc_to_fo = thr_gc_to_fo,
    threshold_fo_to_gc = thr_fo_to_gc,
    dtmax = T_gyro_apex,
    dtmin = 1.0e-2 * T_gyro_0,
    safety_fo = 0.04,
    check_dt = 0.25 * T_gyro_apex,
    check_interval = 4,
    abstol = 1.0e-7,
    reltol = 1.0e-7,
    maxiters = 500_000,
)

alg_curv = AdaptiveHybrid(; mode_kw..., adiabaticity = :curvature)
alg_gradB = AdaptiveHybrid(; mode_kw..., adiabaticity = :gradB)
alg_both = AdaptiveHybrid(; mode_kw..., adiabaticity = :both)
alg_jac = AdaptiveHybrid(; mode_kw..., adiabaticity = :jacobian)

sol_curv = solve(prob_fo, alg_curv; verbose = false, seed = 42).u[1]
sol_gradB = solve(prob_fo, alg_gradB; verbose = false, seed = 42).u[1]
sol_both = solve(prob_fo, alg_both; verbose = false, seed = 42).u[1]
sol_jac = solve(prob_fo, alg_jac; verbose = false, seed = 42).u[1];

function fo_time_frac(s)
    ad = s.stats.adiabaticity
    isempty(ad.t) && return 0.0
    t_end = s.t[end]
    t_fo = 0.0
    for i in eachindex(ad.t)
        t_hi = i < length(ad.t) ? ad.t[i + 1] : t_end
        ad.mode[i] === :FO && (t_fo += t_hi - ad.t[i])
    end
    return t_fo / (t_end - ad.t[1])
end

function rms_diff(s)
    errs = [norm(u[1:3] - sol_fo(t)[1:3]) for (t, u) in zip(s.t, s.u)]
    return sqrt(mean(errs .^ 2))
end

io2 = IOBuffer() #hide
println(io2, "| Mode | FO time fraction | Saved pts | RMS pos. diff [m] |") #hide
println(io2, "| :--- | ---: | ---: | ---: |") #hide
for (name, s) in ( #hide
        ("`:curvature`", sol_curv), #hide
        ("`:gradB`", sol_gradB), #hide
        ("`:both`", sol_both), #hide
        ("`:jacobian`", sol_jac), #hide
    ) #hide
    Printf.@printf( #hide
        io2, "| %s | %.2f | %d | %.2f |\n", #hide
        name, fo_time_frac(s), length(s.t), rms_diff(s), #hide
    ) #hide
end #hide
Markdown.parse(String(take!(io2))) #hide

# Plotting the stored adiabaticity trace for each criterion confirms how `:gradB` and
# `:jacobian` remain stricter in the converging strong-field legs while all four criteria
# trigger full-orbit integration across every apex crossing:

f_comp = Figure(; size = (980, 440), fontsize = 20)
ax_comp = Axis(
    f_comp[1, 1];
    xlabel = L"t / T_{\mathrm{gyro,apex}}",
    ylabel = L"\varepsilon",
    yscale = log10,
    title = "Selected adiabaticity criterion along the curved loop",
    limits = ((t_norm[1], t_norm[end]), (ε_clamp_lo, 1.0)),
)
for (lbl, s, col, ls) in (
        (":curvature", sol_curv, wc[1], :solid),
        (":gradB", sol_gradB, wc[2], :solid),
        (":both", sol_both, :black, :dash),
        (":jacobian", sol_jac, wc[4], :solid),
    )
    ad = s.stats.adiabaticity
    lines!(
        ax_comp, ad.t ./ T_gyro_apex, clamp.(ad.components, ε_clamp_lo, 1.0);
        color = col, linestyle = ls, linewidth = 2.0, label = lbl,
    )
end
hlines!(
    ax_comp, [thr_gc_to_fo];
    color = :gray40, linestyle = :dot, linewidth = 1.6, label = "α = 0.10",
)
axislegend(ax_comp; position = :rb, nbanks = 3, backgroundcolor = (:white, 0.9))

f_comp = DisplayAs.PNG(f_comp) #hide
