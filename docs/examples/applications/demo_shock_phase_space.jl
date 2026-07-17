# # Shock Phase Space
#
# This example demonstrates how to trace ions across a collisionless shock and analyze their phase space distribution, inspired by the demo from IRF-matlab.
# We utilize Liouville's theorem (phase space density conservation), backward/forward tracing, and flux injection to reconstruct the distribution function.

import DisplayAs #hide
using TestParticle
import TestParticle as TP
using StaticArrays
using Random
using FHist
using VelocityDistributionFunctions
using CairoMakie
using Meshes
CairoMakie.activate!(type = "png") #hide

seed = 42;

# ## Upstream Plasma Parameters

const T_ion = 20.0  # ion temperature [eV]
const vth_ion = sqrt(2 * TP.qᵢ * T_ion / TP.mᵢ) # ion thermal speed [m/s]
const V_sw = -400.0e3 # solar wind bulk speed [m/s]
const P_sw = 0.08e-9; # solar wind dynamic pressure [Pa]

# ## Shock Structure Parameters

const n_up = 3.0e6 # upstream number density [m⁻³]
const n_down = 8.0e6 # downstream number density [m⁻³]
const shock_width = 5.0e3; # shock ramp width [m]

# ## Magnetic Field Parameters

const θ_Bn = 45.0 # shock normal angle [degree]
const B_mag = 30.0e-9 # upstream magnetic field magnitude [T]

function compute_tanh_profile_coefficients(θ_Bn, B_mag)
    B_up_y, B_up_x = B_mag .* sincosd(θ_Bn)
    B_up_mag = B_mag

    B_down_x = B_up_x
    B_down_y = 3 * B_up_y
    B_down_mag = sqrt(B_down_x^2 + B_down_y^2)

    B_jump = 0.5 * (B_down_mag - B_up_mag)
    B_avg = 0.5 * (B_up_mag + B_down_mag)
    return B_jump, B_avg
end

const B_jump, B_avg = compute_tanh_profile_coefficients(θ_Bn, B_mag)
const B_normal = 5.0e-9; # shock normal component of B [T]

# ## Field Definitions
# We define custom analytical functions for the electric and magnetic fields across the shock transition layer.

function get_B_shock(r)
    x = r[1]
    bx = B_normal
    by = -B_jump * tanh(x / shock_width) + B_avg
    bz = 0.0
    return SVector{3}(bx, by, bz)
end

"""
Electric field from generalized Ohm's law (Hall term + electron pressure).
"""
function get_E_shock(r)
    xnorm = r[1] / shock_width
    tanh_v = tanh(xnorm)
    sech_v = sech(xnorm)

    ni = -n_up * tanh_v + n_down
    jz = -B_jump * sech_v^2 / (TP.μ₀ * shock_width) # Ampere's law

    by = -B_jump * tanh_v + B_avg
    eni = TP.qᵢ * ni

    ex = -jz * by / eni + P_sw * sech_v^2 / (eni * shock_width)
    ey = jz * B_normal / eni
    ez = -V_sw * (B_avg - B_jump)

    return SVector{3}(ex, ey, ez)
end;

# ## Simulation Setup

nparticles = 10000
const x_source = SA[300.0e3, 0.0, 0.0] # source plane location [m]
const tspan = (0.0, 20.0) # simulation time span [s]
const dt = get_gyroperiod(3 * B_mag) / 20 # time step [s]

param = prepare(get_E_shock, get_B_shock; species = Proton)

## Source velocity distribution (isotropic Maxwellian)
const p_thermal = n_up * TP.qᵢ * T_ion
const vdf = TP.Maxwellian(SA[V_sw, 0.0, 0.0], p_thermal, n_up; m = TP.mᵢ)

function prob_func_maxwellian(prob, ctx)
    v = rand(ctx.rng, vdf)
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

u0_dummy = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func_maxwellian)

println("Starting simulation with $nparticles particles...")
t_mc = @elapsed sols = TP.solve(
    prob, Boris(); dt, savestepinterval = 10, trajectories = nparticles, seed
);
println("Simulation complete. Flux injection tracing time: $(round(t_mc; digits = 2)) s")

## Detector planes (upstream and downstream of the shock)
const x_upstream = 2.0e5  # [m]
const x_downstream = -2.0e5 # [m]
detector_up = Meshes.Plane(
    Meshes.Point(x_upstream, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0)
)
detector_down = Meshes.Plane(
    Meshes.Point(x_downstream, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0)
);

# To get the velocity space distributions, we bin the crossing events into 2D orthogonal velocity planes, integrating over the third dimension.
#
# ## Inputs and Outputs of the Three Methods
# Although the three methods produce similar-looking velocity-space projections, they differ
# substantially in what they take as input and what they return:
#
# | Method | Input (at the **source** unless noted) | Output (at the detector plane) |
# | :--- | :--- | :--- |
# | **1. Forward Monte-Carlo Injection** | Macro-particles launched from `x_source` with initial velocities sampled from the source `Maxwellian` (`vdf`). Each crossing carries its **initial** `v_x` for kinematic flux weighting. | 2-D projected phase-space **density** (flux-weighted histograms) at the upstream/downstream detector, in `[s²/km⁵]`. One value per detector crossing. |
# | **2. Forward Liouville Tracking** | A uniform **sphere** of initial velocities at `x_source` (radius `3 vth_ion`); each sample is weighted by `n0·pdf(vdf)` evaluated at the source (analytical weight). | 2-D projected phase-space **density** (Liouville-weighted histograms) at the upstream/downstream detector, in `[s³/km⁶]` before projection (`[s²/km⁵]` after). |
# | **3. Backward Liouville Tracing** | A regular **velocity grid** placed at the **detector** plane (`vx, vy, vz`); each grid point is traced *backward* in time to `x_source`. | 3-D phase-space **density** `f_3d` evaluated at the detector as `n0·pdf(vdf, traced-source-state)` (gridded, in `[s³/km⁶]`); 2-D projections are obtained by summing over the third axis. |
#
# Key contrasts:
# - Methods 1 and 2 are **source-sampled**: they start at `x_source` and integrate *forward*.
# - Method 3 is **target-sampled**: it starts at the detector and integrates *backward*, so the
#   detector distribution is filled uniformly over the chosen grid regardless of how few source
#   particles would have landed there.
# - Methods 1 and 2 return binned histograms (statistical noise ∝ 1/√N); Method 3 returns a
#   deterministic grid evaluation (no statistical noise, only the grid resolution limit).
#
# ## Method 1: Forward Monte-Carlo Injection
# Simulated particles are treated as macro-particles launched in a steady-state flux at the source.
# Crossing events are converted to phase-space density via kinematic weighting (see below).

function reconstruct_flux_projections(sols, detector, n0, dv_km)
    ## Initial velocities at the source plane
    vxi = [s.u[1][4] for s in sols.u] # initial vx [m/s]
    ## Detect crossings at the plane
    vs, ws_init = get_particle_crossings(sols, detector, vxi)

    v_edges = -1000:dv_km:1000
    h_3d = Hist3D(; binedges = (v_edges, v_edges, v_edges))

    ## Flux normalization factor S = n0_km3 / (N_total * dv_km^2)
    ## Conversion: n0 [m^-3] * 1e9 = n0 [km^-3]
    S = (n0 * 1.0e9) / (length(sols.u) * dv_km^2)

    for (v, vxi_val) in zip(vs, ws_init)
        ## Weight w = (v_xi / v_det) * S
        w = abs(vxi_val) / abs(v[1]) * S
        push!(h_3d, v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3, w) # units: [s^2/km^5]
    end

    return project(h_3d, :z), project(h_3d, :y), project(h_3d, :x)
end

function plot_shock_vdf(hists_up, hists_down, x_up, x_down; vlim = 1000.0)
    fig = Figure(size = (1300, 650), fontsize = 22)
    xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
    ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]

    for i in 1:3
        ## Upstream (row 1) and Downstream (row 2)
        for (row, hists, label, xloc) in
            [(1, hists_up, "Upstream", x_up), (2, hists_down, "Downstream", x_down)]
            ax = Axis(
                fig[row, i], title = "$(label) x = $(xloc * 1.0e-3) km",
                xlabel = xlabels[i], ylabel = ylabels[i];
                xlabelsize = 26, ylabelsize = 26, titlesize = 24,
                xticklabelsize = 20, yticklabelsize = 20,
                limits = (-vlim, vlim, -vlim, vlim)
            )
            h = hists[i]
            hm = h isa Tuple ? heatmap!(ax, h...; colormap = :turbo) :
                heatmap!(ax, h; colormap = :turbo)
            if i == 3
                Colorbar(
                    fig[row, 4], hm; label = L"[\mathrm{s}^2/\mathrm{km}^5]",
                    labelsize = 22, ticklabelsize = 18
                )
            end
        end
    end
    return fig
end

function plot_downstream_comparison(h1, h2, h3; vlim = 1000.0)
    titles = ["Flux Injection", "Forward Liouville", "Backward Liouville"]
    cols = [L"V_x–V_y", L"V_x–V_z", L"V_y–V_z"]
    xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
    ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]
    hists = (h1, h2, h3)
    fig = Figure(size = (1300, 1100), fontsize = 22)
    gl = fig[1, 1] = GridLayout()
    Label(gl[1, 2:4], "Downstream velocity distributions (x = -200 km)"; fontsize = 24, tellwidth = false)
    for r in 1:3
        for i in 1:3
            ax = Axis(
                gl[r + 1, i + 1], title = cols[i],
                xlabel = xlabels[i], ylabel = ylabels[i];
                xlabelsize = 22, ylabelsize = 22, titlesize = 20,
                xticklabelsize = 18, yticklabelsize = 18,
                limits = (-vlim, vlim, -vlim, vlim),
                aspect = 1
            )
            h = hists[r][i]
            hm = h isa Tuple ? heatmap!(ax, h...; colormap = :turbo) :
                heatmap!(ax, h; colormap = :turbo)
            if i == 3
                Colorbar(gl[r + 1, 5], hm; labelsize = 18, ticklabelsize = 14)
            end
        end
        Label(gl[r + 1, 1], titles[r]; fontsize = 18, rotation = π / 2, tellheight = false)
    end
    return fig
end

hists_up = reconstruct_flux_projections(sols, detector_up, n_up, 20.0)
hists_down = reconstruct_flux_projections(sols, detector_down, n_up, 20.0)

fig_flux = plot_shock_vdf(hists_up, hists_down, x_upstream, x_downstream)
fig_flux = DisplayAs.PNG(fig_flux) #hide

# The kinematic weight ``w = |v_{x,\mathrm{init}}| / |v_{x,\mathrm{det}}|`` converts the density-sampled
# launch into a steady-state flux (``|v_{x,\mathrm{init}}|``) and back to phase-space density at the
# detector (``1/|v_{x,\mathrm{det}}|``). Dropping both factors (equal-weight binning) leaves a small
# upstream error, where particles are nearly unperturbed, but downstream it noticeably underestimates
# the density in the low-``|v_x|`` tails of reflected and decelerated ions.
#
# ## Method 2: Forward Liouville Tracking
# Forward Liouville tracking starts from a sphere of initial conditions in velocity space at the source
# and traces forward to the detector, combining Monte-Carlo sampling of the sphere with Liouville's theorem.

function reconstruct_liouville_projections(sols, detector, vdf, n0, Vsphere; dv_km = 20.0)
    ## 1. Initial weights from source PDF
    ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols.u]
    ## 2. Crossings
    vs, ws = get_particle_crossings(sols, detector, ws0)

    v_edges = -1000:dv_km:1000
    h_3d = Hist3D(; binedges = (v_edges, v_edges, v_edges))

    ## Normalization in km-based units: Vsphere [km^3], dv_km [km/s]
    ## Conversion: 1 m^3 = 1e-9 km^3
    Vsphere_km = Vsphere * 1.0e-9
    S_L = Vsphere_km / (length(sols.u) * dv_km^2)

    for (v, w) in zip(vs, ws)
        ## w is in [s^3/m^6]. Convert to [s^3/km^6] by multiplying 1e18.
        push!(h_3d, v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3, w * 1.0e18 * S_L)
    end

    return project(h_3d, :z), project(h_3d, :y), project(h_3d, :x)
end

nparticles_m2 = 10000
const vradius_m2 = 3 * vth_ion # velocity space radius, [m/s]
const Vsphere_m2 = (4 / 3) * π * vradius_m2^3 # velocity space volume

## Uniform sampling in a 3D sphere
function prob_func_m2(prob, ctx)
    r = vradius_m2 * rand(ctx.rng)^(1 / 3)
    ϕ = 2π * rand(ctx.rng)
    θ = acos(2 * rand(ctx.rng) - 1)

    sinθ, cosθ = sincos(θ)
    cosϕ, sinϕ = sincos(ϕ)
    v = SA[V_sw + r * sinθ * cosϕ, r * sinθ * sinϕ, r * cosθ]
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

prob_m2 = TraceProblem(
    SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param; prob_func = prob_func_m2
)
t_liou = @elapsed sols_m2 = TP.solve(
    prob_m2, Boris(); dt, savestepinterval = 10, trajectories = nparticles_m2, seed
);

hists_up_m2 = reconstruct_liouville_projections(
    sols_m2, detector_up, vdf, n_up, Vsphere_m2
)
hists_down_m2 = reconstruct_liouville_projections(
    sols_m2, detector_down, vdf, n_up, Vsphere_m2
)

fig_forward = plot_shock_vdf(hists_up_m2, hists_down_m2, x_upstream, x_downstream)
fig_forward = DisplayAs.PNG(fig_forward) #hide

# ## Method 3: Backward Tracing
# Starting from a velocity-space grid at the detector, particles are traced backward; the phase-space
# density at the detector equals the source density evaluated at the traced initial state. Integrating
# the 3D values gives the projections compared against the other methods below.
#
# The grid uses `dv_km = 20` km/s, with `vy_range = ±400` km/s to cover the band around `V_sw = -400`
# km/s and `vx`/`vz = ±1000` km/s to capture the reflected beam at positive velocities. Since the
# populated region is a small fraction of the box, `adaptive = true` first traces a coarse grid
# (`dv_coarse_km = 60` km/s) to locate the active cells, then a fine grid only inside that box (padded
# by `margin_km`). Discarded cells lie below the `f_max·10⁻⁶` clipping threshold, so the projections
# are unchanged while the trajectory count drops by ≈2–3×.
#
# **Termination note.** A back-traced trajectory is stopped as soon as it (i) moves *away* from the
# source (`u[1] < detector_x - 100 km`) or (ii) has passed the source plane by a safe margin
# (`u[1] > x_source[1] + 6000 km`). The margin must exceed the distance covered in one saved step at the
# highest grid speed (≈5.3×10³ km here), otherwise the post-crossing saved point is missed and the cell
# is dropped. A too-tight guard beyond the source formerly truncated the trajectory to a single point and
# removed the entire fast half of the upstream Maxwellian in the `V_x–V_z` plot ("cut in half"). Stopping
# at the source also bounds the integration length (≈half the steps) without changing the recovered
# density.

## Solve one backward-tracing pass over a uniform velocity grid and return the 3D phase-space
## density sampled at the source plane, in [s³/km⁶].
function run_backward_pass(vx_grid, vy_grid, vz_grid, detector_x, vdf, n0, dt, param)
    nx, ny, nz = length(vx_grid), length(vy_grid), length(vz_grid)
    ntraj = nx * ny * nz

    function prob_func(prob, ctx)
        iz = (ctx.sim_id - 1) % nz + 1
        iy = ((ctx.sim_id - 1) ÷ nz) % ny + 1
        ix = ((ctx.sim_id - 1) ÷ (nz * ny)) % nx + 1
        u0 = SA[detector_x, 0.0, 0.0, vx_grid[ix], vy_grid[iy], vz_grid[iz]]
        return remake(prob, u0 = u0)
    end

    source_plane = Meshes.Plane(Meshes.Point(x_source...), Meshes.Vec(1.0, 0.0, 0.0))
    prob = TraceProblem(
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], (0.0, -20.0), param;
        prob_func = prob_func
    )

    sols = TP.solve(
        prob, Boris(), EnsembleThreads(); dt = -dt, trajectories = ntraj,
        savestepinterval = 10,
        ## Stop at the source plane (and terminate away-going trajectories). See "Method 3".
        isoutside = (u, p, t) -> u[1] < detector_x - 1.0e5 || u[1] > x_source[1] + 6000.0e3
    )

    f_3d = zeros(nx, ny, nz)
    for (i, sol) in enumerate(sols.u)
        st = get_first_crossing(sol, source_plane)
        if !any(isnan, st)
            iz = (i - 1) % nz + 1
            iy = ((i - 1) ÷ nz) % ny + 1
            ix = ((i - 1) ÷ (nz * ny)) % nx + 1
            f_3d[ix, iy, iz] = n0 * pdf(vdf, st[SA[4, 5, 6]]) * 1.0e18
        end
    end
    return f_3d
end

function reconstruct_backward_projections(
        detector_x, vdf, n0, dt, param;
        v_range = 1000.0e3, vy_range = 400.0e3, dv_km = 20.0,
        adaptive = true, dv_coarse_km = 60.0, margin_km = 150.0
    )
    dv = dv_km * 1.0e3

    if adaptive
        ## Pass 1 (coarse): locate the populated region cheaply.
        vx_c = range(-v_range, v_range, step = dv_coarse_km * 1.0e3)
        vy_c = range(-vy_range, vy_range, step = dv_coarse_km * 1.0e3)
        vz_c = range(-v_range, v_range, step = dv_coarse_km * 1.0e3)
    else
        ## Uniform grid over the full box (baseline / no adaptation).
        vx_c = range(-v_range, v_range, step = dv)
        vy_c = range(-vy_range, vy_range, step = dv)
        vz_c = range(-v_range, v_range, step = dv)
    end

    t_solve = @elapsed begin
        f_coarse = run_backward_pass(vx_c, vy_c, vz_c, detector_x, vdf, n0, dt, param)
        if adaptive
            ## Keep every cell well above the display threshold, padded by `margin_km`.
            kept = findall(f_coarse .> maximum(f_coarse) * 1.0e-5)
            if isempty(kept)
                vx_grid, vy_grid, vz_grid = vx_c, vy_c, vz_c
                f_3d_km = f_coarse
            else
                ixs, iys, izs = getindex.(kept, 1), getindex.(kept, 2), getindex.(kept, 3)
                vx_grid = range(
                    max(-v_range, floor((vx_c[minimum(ixs)] - margin_km * 1.0e3) / dv) * dv),
                    min(v_range, ceil((vx_c[maximum(ixs)] + margin_km * 1.0e3) / dv) * dv);
                    step = dv
                )
                vy_grid = range(
                    max(-vy_range, floor((vy_c[minimum(iys)] - margin_km * 1.0e3) / dv) * dv),
                    min(vy_range, ceil((vy_c[maximum(iys)] + margin_km * 1.0e3) / dv) * dv);
                    step = dv
                )
                vz_grid = range(
                    max(-v_range, floor((vz_c[minimum(izs)] - margin_km * 1.0e3) / dv) * dv),
                    min(v_range, ceil((vz_c[maximum(izs)] + margin_km * 1.0e3) / dv) * dv);
                    step = dv
                )
                f_3d_km = run_backward_pass(vx_grid, vy_grid, vz_grid, detector_x, vdf, n0, dt, param)
            end
        else
            vx_grid, vy_grid, vz_grid = vx_c, vy_c, vz_c
            f_3d_km = f_coarse
        end
    end
    nparticles_bw = length(vx_grid) * length(vy_grid) * length(vz_grid)

    ## Integral over third dimension [km/s]: f_int in [s²/km⁵]
    f_xy = dropdims(sum(f_3d_km, dims = 3), dims = 3) .* (step(vz_grid) * 1.0e-3)
    f_xz = dropdims(sum(f_3d_km, dims = 2), dims = 2) .* (step(vy_grid) * 1.0e-3)
    f_yz = dropdims(sum(f_3d_km, dims = 1), dims = 1) .* (step(vx_grid) * 1.0e-3)

    for f in (f_xy, f_xz, f_yz)
        f_max = maximum(f)
        for i in eachindex(f)
            if f[i] < f_max * 1.0e-6
                f[i] = NaN
            end
        end
    end

    return (
            (vx_grid .* 1.0e-3, vy_grid .* 1.0e-3, f_xy),
            (vx_grid .* 1.0e-3, vz_grid .* 1.0e-3, f_xz),
            (vy_grid .* 1.0e-3, vz_grid .* 1.0e-3, f_yz),
        ), t_solve, nparticles_bw
end

res_up_bw, t_bw_up, n_bw_up =
    reconstruct_backward_projections(x_upstream, vdf, n_up, dt, param)
res_down_bw, t_bw_down, n_bw_down =
    reconstruct_backward_projections(x_downstream, vdf, n_up, dt, param)
t_bw = t_bw_up + t_bw_down
n_bw = n_bw_up + n_bw_down

fig_backward = plot_shock_vdf(res_up_bw, res_down_bw, x_upstream, x_downstream)
fig_backward = DisplayAs.PNG(fig_backward) #hide

# The reconstructed phase-space density on both sides of the shock.
#
# We then compare the downstream distribution with the two forward estimates.
fig_cmp = plot_downstream_comparison(hists_down, hists_down_m2, res_down_bw)
fig_cmp = DisplayAs.PNG(fig_cmp) #hide

# ## Summary
# This example illustrates three complementary ways to reconstruct the phase space density from particle simulations.
#
# ### Cost and which method to use
# The wall-clock times in the table below are an **apples-to-oranges** comparison: Methods 1–2 launch
# 10⁴ macro-particles to build a *noisy* histogram, whereas Method 3 fills a *complete, deterministic*
# velocity grid of ≈3×10⁵ cells. To reproduce the same grid with forward tracing you would need roughly
# (cells) × (≈10²–10³ particles per cell for usable statistics) ≈ 10⁷–10⁹ trajectories — orders of
# magnitude slower. So backward Liouville *is* the most **cost-effective** method once you require a full,
# noise-free `f(v)`: it is target-sampled (no particles wasted in empty phase-space) and deterministic (no
# Monte-Carlo oversampling).
#
# Per trajectory it is also a few × cheaper (shorter effective integration, no statistical weighting, and
# away-going trajectories terminate at once). Stopping each back-traced trajectory right after the source
# plane (the `isoutside` guard in Method 3) cuts the integration length roughly in half; in this `dt`
# regime the wall time is dominated by per-trajectory solve/threading overhead, so the saving is modest,
# but the bound is worthwhile and changes no result.
#
# **Guidance.**
# - **Method 1 (flux injection):** quick statistical picture or low-order moments / flux; cheapest but noisy, with poor tail resolution.
# - **Method 2 (forward Liouville):** better than MC without a full grid; still source-sampled, so empty regions are wasted.
# - **Method 3 (backward Liouville):** a complete, noise-free distribution on a *chosen* grid (best accuracy per cell); pair with the adaptive two-pass and source-plane termination to keep the trajectory count down. This is the method to reach for when accuracy matters.

t_per_mc = t_mc / nparticles * 1.0e6
t_per_liou = t_liou / nparticles_m2 * 1.0e6
t_per_bw = t_bw / n_bw * 1.0e6;

using Markdown, Printf #hide
io = IOBuffer() #hide
println(io, "| Method | Flux Injection | Forward Liouville Tracking | Backward Liouville Tracing |") #hide
println(io, "| :--- | :--- | :--- | :--- |") #hide
println(io, "| **Noise** | Statistical (∝ 1/√N) | Low (analytical weights) | None (grid-based) |") #hide
println(io, "| **Coverage** | Source-sampled | Source-sampled | Target-sampled |") #hide
println(io, "| **Tail resolution** | Poor without large N | Limited by sphere radius | Uniform across grid |") #hide
println(io, "| **Post-processing** | Binning + weighting | Binning + projection | PDF evaluation only |") #hide
@printf( #hide
    io, "| **Cost** | %.1f s (%.1f µs/traj, %d traj.) | %.1f s (%.1f µs/traj, %d traj.) | %.1f s (%.1f µs/traj, %d traj.) |\n", #hide
    t_mc, t_per_mc, nparticles, t_liou, t_per_liou, nparticles_m2, t_bw, t_per_bw, n_bw #hide
) #hide
Markdown.parse(String(take!(io))) #hide
