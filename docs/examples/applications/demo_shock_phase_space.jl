# # Shock Phase Space
#
# This example demonstrates how to trace ions across a collisionless shock and analyze their
# phase space distribution, inspired by the demo from IRF-matlab.
# We utilize Liouville's theorem (phase space density conservation), backward/forward tracing,
# and flux injection to reconstruct the distribution function.

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
# Custom analytical electric and magnetic fields across the shock transition layer.

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
#
# The source plane is placed just upstream of the shock, far enough that the fastest particle
# gyroradius (~350 km at 1000 km/s in 30 nT) fits inside the uniform upstream region without
# the source plane clipping any gyro-orbit.

nparticles = 10000
const x_source = SA[500.0e3, 0.0, 0.0] # source plane location [m]
const tspan = (0.0, 20.0) # forward simulation time span [s]
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
    prob, Boris(); dt, savestepinterval = 1, trajectories = nparticles, seed
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

# To get the velocity space distributions, we bin the crossing events into 2D orthogonal
# velocity planes, integrating over the third dimension.
#
# ## What each method computes
#
# All three methods return the **phase-space density** ``f`` at the detector, in
# ``[\mathrm{s}^3/\mathrm{km}^6]`` (3D) or ``[\mathrm{s}^2/\mathrm{km}^5]`` (2D projection):
#
# | Method | Input | Output |
# | :--- | :--- | :--- |
# | **1. Forward Monte-Carlo** | Macro-particles launched from `x_source` with velocities sampled from the source `Maxwellian` (`vdf`). Each crossing contributes a constant weight `S = n0_km³ / (N · dv²)`, so the histogram directly estimates ``f``. | 2-D projected ``f`` (histogram), `[s²/km⁵]`. |
# | **2. Forward Liouville** | A uniform **sphere** of initial velocities at `x_source`; each sample carries weight `n0·pdf(vdf, v_source)` (the source ``f``). By Liouville's theorem `f_det(v_det) = f_source(v_source)`, so binning detector velocities with the source ``f`` value gives ``f`` at the detector. | 2-D projected ``f`` (bin-averaged), `[s²/km⁵]`. |
# | **3. Backward Liouville** | A regular **velocity grid** at the **detector**; each grid point is traced *backward* to `x_source` and `f_det = n0·pdf(vdf, v_traced)` is evaluated. | 3-D ``f`` on a grid `[s³/km⁶]`; 2-D projections by summing. |

# ## Method 1: Forward Monte-Carlo
#
# Particles are launched from the source with velocities drawn from the source Maxwellian.
# Each detector crossing contributes a constant weight `S`, so the histogram estimates the
# phase-space density ``f`` at the detector.

function reconstruct_flux_projections(sols, detector, n0, dv_km)
    vs, _ = get_particle_crossings(sols, detector, 1.0)

    v_edges = -1000:dv_km:1000
    h_3d = Hist3D(; binedges = (v_edges, v_edges, v_edges))

    ## S = n0 [km⁻³] / (N_total · dv_km²)  →  [s²/km⁵] after 1D projection
    S = (n0 * 1.0e9) / (length(sols.u) * dv_km^2)

    for v in vs
        push!(h_3d, v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3, S)
    end

    return project(h_3d, :z), project(h_3d, :y), project(h_3d, :x)
end

# ## Plotting helpers
#
# A single log-scale colour map is used everywhere.  `compute_common_cr` finds a global
# `(fmin, fmax)` across a set of 2-D distributions so that the comparison plot uses one
# shared colour range, making physical differences visible rather than hidden by per-panel
# scaling.

"""
Display a 2-D histogram or (x, y, A) tuple on a log10 colour scale.
Zeros / NaNs are clamped to `fmin` so empty cells use the bottom of the colour bar.
"""
function logheatmap!(ax, h::Union{Hist2D, Tuple}; colormap = :turbo, cr = (1e-8, 1.0))
    if h isa Tuple
        x, y, A = h
    else
        edges = collect.(binedges(h))
        x = 0.5 .* (edges[1][2:end] .+ edges[1][1:end-1])
        y = 0.5 .* (edges[2][2:end] .+ edges[2][1:end-1])
        A = h.bincounts
    end
    fmin, fmax = cr
    B = float(A)
    replace!(x -> isfinite(x) && x > 0 ? x : fmin, B)
    hm = heatmap!(ax, x, y, B; colormap, colorscale = log10, colorrange = (fmin, fmax))
    return hm
end

"""
Return (fmin, fmax) for a shared log colour range across several 2-D distributions.
`rel` sets the floor as a fraction of the global maximum.
"""
function compute_common_cr(hists; rel = 1.0e-5)
    fmax = 0.0
    for h in hists
        A = h isa Tuple ? h[3] : h.bincounts
        m = maximum(A)
        m > fmax && (fmax = m)
    end
    fmax ≤ 0 && return (1.0e-30, 1.0)
    return (fmax * rel, fmax)
end

function plot_shock_vdf(hists_up, hists_down, x_up, x_down; vlim = 1000.0)
    cr = compute_common_cr((hists_up..., hists_down...))
    fig = Figure(size = (1300, 650), fontsize = 22)
    xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
    ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]

    for i in 1:3
        for (row, hists, label, xloc) in
            [(1, hists_up, "Upstream", x_up), (2, hists_down, "Downstream", x_down)]
            ax = Axis(
                fig[row, i], title = "$(label) x = $(xloc * 1.0e-3) km",
                xlabel = xlabels[i], ylabel = ylabels[i];
                xlabelsize = 26, ylabelsize = 26, titlesize = 24,
                xticklabelsize = 20, yticklabelsize = 20,
                limits = (-vlim, vlim, -vlim, vlim)
            )
            hm = logheatmap!(ax, hists[i]; cr)
            if i == 3
                Colorbar(
                    fig[row, 4], hm; label = L"\log_{10}([\mathrm{s}^2/\mathrm{km}^5])",
                    labelsize = 22, ticklabelsize = 18
                )
            end
        end
    end
    return fig
end

"""
Three-row comparison of the downstream `f` from each method, with a single shared colour bar.
"""
function plot_downstream_comparison(h1, h2, h3; vlim = 1000.0)
    titles = ["Flux Injection", "Forward Liouville", "Backward Liouville"]
    xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
    ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]
    hists = (h1, h2, h3)
    all_h = (h1..., h2..., h3...)
    cr = compute_common_cr(all_h)
    fig = Figure(size = (1300, 1100), fontsize = 22)
    gl = fig[1, 1] = GridLayout()
    Label(
        gl[1, 2:4], "Downstream velocity distributions (x = $(x_downstream*1.0e-3) km)";
        fontsize = 30, tellwidth = false
    )
    for r in 1:3
        row = r + 1
        for i in 1:3
            ax = Axis(
                gl[row, i + 1],
                xlabel = xlabels[i], ylabel = ylabels[i];
                xlabelsize = 28, ylabelsize = 28,
                xticklabelsize = 18, yticklabelsize = 18,
                limits = (-vlim, vlim, -vlim, vlim),
                aspect = 1
            )
            hm = logheatmap!(ax, hists[r][i]; cr)
            if r == 1 && i == 3
                Colorbar(gl[2:4, 5], hm;
                    label = L"\log_{10}([\mathrm{s}^2/\mathrm{km}^5])",
                    labelsize = 18, ticklabelsize = 14)
            end
        end
        Label(gl[row, 1], titles[r]; fontsize = 24, rotation = π / 2, tellheight = false)
    end
    return fig
end

hists_up = reconstruct_flux_projections(sols, detector_up, n_up, 20.0)
hists_down = reconstruct_flux_projections(sols, detector_down, n_up, 20.0)

fig_flux = plot_shock_vdf(hists_up, hists_down, x_upstream, x_downstream)
fig_flux = DisplayAs.PNG(fig_flux) #hide

# Each crossing contributes a constant weight `S = n0 / (N · dv²)`, so the histogram estimates
# the phase-space density ``f`` at the detector (the `N` and `dv` factors are removed).
#
# ## Method 2: Forward Liouville Tracking
#
# Forward Liouville tracking starts from a sphere of initial conditions in velocity space at
# the source and traces forward to the detector.  By Liouville's theorem
# ``f_{\det}(\mathbf{v}_{\det}) = f_{\src}(\mathbf{v}_{\src})``; the source value
# ``n_0\,\mathrm{pdf}(\mathrm{vdf}, \mathbf{v}_{\src})`` is carried unchanged along each
# trajectory.  We bin the *detector* velocities and average the source ``f`` inside each
# bin, exactly mirroring Method 3.

function reconstruct_liouville_projections(sols, detector, vdf, n0; dv_km = 20.0)
    ## Source f for every trajectory [s³/m⁶]
    ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols.u]
    ## Detector crossings carry the source f (Liouville)
    vs, ws = get_particle_crossings(sols, detector, ws0)

    v_edges = -1000:dv_km:1000
    nb = length(v_edges) - 1
    centers = 0.5 .* (v_edges[2:end] .+ v_edges[1:end-1])
    sum_f = zeros(nb, nb, nb)
    count = zeros(Int, nb, nb, nb)

    for (v, w) in zip(vs, ws)
        ix = floor(Int, (v[1] * 1.0e-3 - v_edges[1]) / dv_km) + 1
        iy = floor(Int, (v[2] * 1.0e-3 - v_edges[1]) / dv_km) + 1
        iz = floor(Int, (v[3] * 1.0e-3 - v_edges[1]) / dv_km) + 1
        if 1 <= ix <= nb && 1 <= iy <= nb && 1 <= iz <= nb
            ## [s³/m⁶] → [s³/km⁶]
            sum_f[ix, iy, iz] += w * 1.0e18
            count[ix, iy, iz] += 1
        end
    end
    f_3d = ifelse.(count .> 0, sum_f ./ count, 0.0)

    dv = dv_km * 1.0e3
    f_xy = dropdims(sum(f_3d, dims = 3), dims = 3) .* (dv * 1.0e-3)
    f_xz = dropdims(sum(f_3d, dims = 2), dims = 2) .* (dv * 1.0e-3)
    f_yz = dropdims(sum(f_3d, dims = 1), dims = 1) .* (dv * 1.0e-3)

    return ((centers, centers, f_xy), (centers, centers, f_xz), (centers, centers, f_yz))
end

nparticles_m2 = 100000
const vradius_m2 = 3 * vth_ion # velocity-space radius, [m/s]

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
    prob_m2, Boris(); dt, savestepinterval = 1, trajectories = nparticles_m2, seed
);

hists_up_m2 = reconstruct_liouville_projections(sols_m2, detector_up, vdf, n_up)
hists_down_m2 = reconstruct_liouville_projections(sols_m2, detector_down, vdf, n_up)

fig_forward = plot_shock_vdf(hists_up_m2, hists_down_m2, x_upstream, x_downstream)
fig_forward = DisplayAs.PNG(fig_forward) #hide

# The sphere radius `3 vth_ion` covers the bulk of the Maxwellian; `nparticles_m2 = 10⁵`
# gives usable statistics in the populated bins.  The bin-averaged source ``f`` is a direct
# Monte-Carlo estimate of ``f_{\det}``, on the same grid and with the same units as Method 3,
# so the two should agree up to sampling noise.
#
# ## Method 3: Backward Liouville Tracing
#
# Starting from a velocity-space grid at the detector, each grid point is traced *backward*
# in time to the source plane.  The phase-space density at the detector equals the source
# density evaluated at the traced-back state: ``f_{\det}(\mathbf{v}_{\det}) = n_0\,
# \mathrm{pdf}(\mathrm{vdf}, \mathbf{v}_{\src})``.
#
# Every step is saved (`savestepinterval = 1`) so that no source-plane crossing is missed,
# and a trajectory is terminated only once it has crossed the source plane and moved a safe
# distance beyond it (`u[1] > x_source + margin`).  Gyrating trajectories that temporarily
# move away are not terminated, so every grid cell whose backward trajectory eventually
# crosses the source receives a value.

## One backward-tracing pass over a uniform velocity grid; returns the 3-D phase-space
## density at the detector in [s³/km⁶].
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
    ## 20 s backward span: the slowest relevant crossing completes well within this window,
    ## giving the same result at ~3x lower cost. The post-source margin must exceed the
    ## distance covered in one saved step at the highest grid speed (~1000 km/s · dt ≈ 2.6 km);
    ## 500 km is more than enough to capture the post-crossing point.
    tspan_bw = (0.0, -20.0)
    post_source_margin = 500.0e3
    prob = TraceProblem(
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan_bw, param;
        prob_func = prob_func
    )

    sols = TP.solve(
        prob, Boris(), EnsembleThreads(); dt = -dt, trajectories = ntraj,
        savestepinterval = 1,
        ## Terminate only after the particle has crossed the source plane and moved a
        ## safe distance beyond it, so the crossing point is always saved.  No away-going
        ## guard: gyrating trajectories must be allowed to come back.
        isoutside = (u, p, t) -> u[1] > x_source[1] + post_source_margin
    )

    f_3d = zeros(nx, ny, nz)
    for (i, sol) in enumerate(sols.u)
        st = get_first_crossing(sol, source_plane)
        if !any(isnan, st)
            iz = (i - 1) % nz + 1
            iy = ((i - 1) ÷ nz) % ny + 1
            ix = ((i - 1) ÷ (nz * ny)) % nx + 1
            ## [s³/m⁶] → [s³/km⁶]
            f_3d[ix, iy, iz] = n0 * pdf(vdf, st[SA[4, 5, 6]]) * 1.0e18
        end
    end
    return f_3d
end

function reconstruct_backward_projections(
        detector_x, vdf, n0, dt, param;
        v_range = 1000.0e3, vy_range = 400.0e3, vz_range = 1000.0e3, dv_km = 20.0,
        adaptive = true, dv_coarse_km = 60.0, margin_km = 150.0
    )
    dv = dv_km * 1.0e3

    if adaptive
        ## Pass 1 (coarse): locate the populated region cheaply.
        vx_c = range(-v_range, v_range, step = dv_coarse_km * 1.0e3)
        vy_c = range(-vy_range, vy_range, step = dv_coarse_km * 1.0e3)
        vz_c = range(-vz_range, vz_range, step = dv_coarse_km * 1.0e3)
    else
        vx_c = range(-v_range, v_range, step = dv)
        vy_c = range(-vy_range, vy_range, step = dv)
        vz_c = range(-vz_range, vz_range, step = dv)
    end

    t_solve = @elapsed begin
        f_coarse = run_backward_pass(vx_c, vy_c, vz_c, detector_x, vdf, n0, dt, param)
        if adaptive
            kept = findall(f_coarse .> maximum(f_coarse) * 1.0e-6)
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
                    max(-vz_range, floor((vz_c[minimum(izs)] - margin_km * 1.0e3) / dv) * dv),
                    min(vz_range, ceil((vz_c[maximum(izs)] + margin_km * 1.0e3) / dv) * dv);
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

# ## Comparison of the three methods (shared colour scale)
#
# The downstream `f` from the three methods is compared using a **single shared colour bar**,
# so the colour ranges are directly comparable.

fig_cmp = plot_downstream_comparison(hists_down, hists_down_m2, res_down_bw)
fig_cmp = DisplayAs.PNG(fig_cmp) #hide

# ## Summary
#
# This example illustrates three complementary ways to reconstruct the phase-space density
# from particle simulations, all returning the same physical quantity ``f`` and sharing a
# common colour scale in the comparison plot.

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
