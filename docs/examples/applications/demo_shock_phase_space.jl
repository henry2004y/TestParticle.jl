# # Shock Phase Space
#
# This example demonstrates how to trace ions across a collisionless shock and analyze their
# phase space distribution, inspired by the demo from IRF-matlab.
# We utilize Liouville's theorem (phase space density conservation), backward/forward tracing,
# and Monte Carlo sampling to reconstruct the distribution function.

import DisplayAs #hide
using TestParticle
import TestParticle as TP
using StaticArrays
using Random
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
## tanh profile from n_up (x → +∞) to n_down (x → −∞)
const n_jump = 0.5 * (n_down - n_up) # half density jump [m⁻³]
const n_avg = 0.5 * (n_down + n_up) # mean density [m⁻³]
const shock_width = 5.0e3; # shock ramp width [m]
## B_t,down / B_t,up, kept below the Rankine-Hugoniot value n_down / n_up ≈ 2.7: at this ramp
## width that jump gives a Hall term above the ion ram energy and reflects the entire beam.
const r_Bt = 2.24

# ## Magnetic Field Parameters

const B_normal = 5.0e-9 # shock normal component of B [T], continuous across the ramp
const B_mag = 30.0e-9 # upstream magnetic field magnitude [T]
## The shock normal is x̂, so θ_Bn follows from the two field strengths above.
const θ_Bn = acosd(B_normal / B_mag) # shock normal angle [degree]
println("Shock normal angle θ_Bn = $(round(θ_Bn; digits = 1))°")

"""
Tanh coefficients of the tangential field: `B_y(x) = -B_jump·tanh(x/w) + B_avg` runs from
`B_mag·sind(θ_Bn)` upstream to `r_Bt` times that value downstream.
"""
function compute_tanh_profile_coefficients(θ_Bn, B_mag, r_Bt)
    B_up_y = B_mag * sind(θ_Bn)
    B_down_y = r_Bt * B_up_y

    B_jump = 0.5 * (B_down_y - B_up_y)
    B_avg = 0.5 * (B_down_y + B_up_y)
    return B_jump, B_avg
end

const B_jump, B_avg = compute_tanh_profile_coefficients(θ_Bn, B_mag, r_Bt)

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

    ni = -n_jump * tanh_v + n_avg
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
## Source phase-space density [s³/m⁶]; `pdf` returns a normalized VDF.
f_src(v) = n_up * pdf(vdf, v)

function prob_func_maxwellian(prob, ctx)
    v = rand(ctx.rng, vdf)
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

u0_dummy = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func_maxwellian)

println("Starting simulation with $nparticles particles...")
t_mc = @elapsed sols = TP.solve(
    prob, Boris(), EnsembleThreads(); dt, savestepinterval = 1,
    trajectories = nparticles, seed
);
println("Simulation complete. Monte Carlo tracing time: $(round(t_mc; digits = 2)) s")

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
# | **1. Forward Monte Carlo** | Macro-particles launched from `x_source` with velocities sampled from the source `Maxwellian` (`vdf`), i.e. density weighted. Each crossing is weighted by `S · \|v_x,src\| / \|v_x,det\|` with `S = n0_km³ / (N · dv²)`: the first factor makes the ensemble flux weighted, the second converts the crossing flux back into a density. | 2-D projected ``f`` (histogram), `[s²/km⁵]`. |
# | **2. Forward Liouville** | A uniform **sphere** of initial velocities at `x_source`; each sample carries the source ``f`` (`n0·pdf(vdf, v_source)`) *and* the velocity volume `vsphere/N` it represents. By Liouville's theorem `f_det(v_det) = f_source(v_source)`, so each crossing deposits `f·ΔV` into the detector bin it lands in, with `ΔV = (vsphere/N)·\|v_x,src\|/\|v_x,det\|`. Summing `f·ΔV` and dividing by the bin volume gives the bin-averaged ``f``. | 2-D projected ``f`` (histogram), `[s²/km⁵]`. |
# | **3. Backward Liouville** | A regular **velocity grid** at the **detector**; each grid point is traced *backward* to `x_source` and `f_det = n0·pdf(vdf, v_traced)` is evaluated. | 3-D ``f`` on a grid `[s³/km⁶]`; 2-D projections by summing. |

# ## Method 1: Forward Monte Carlo
#
# Particles are launched from the source with velocities drawn from the source Maxwellian,
# so the ensemble is density weighted.  A steady beam, however, crosses the source plane
# flux weighted: faster particles are injected more often.  Multiplying each sample by
# ``|v_{x,\mathrm{src}}|`` supplies that weighting, and dividing by ``|v_{x,\det}|`` at the
# detector undoes the flux factor of the recorded crossings.  The net weight
# ``S\,|v_{x,\mathrm{src}}|/|v_{x,\det}|`` reduces to a constant only when every particle
# crosses the detector once with an unchanged ``v_x`` — which is why the simple constant
# weight is exact upstream but not downstream of the shock.

function reconstruct_mc_projections(sols, detector, n0, dv_km)
    ## The launch velocities are drawn from the source VDF, i.e. density weighted, whereas a
    ## steady beam crosses the source plane flux weighted. Carrying |v_x,src| as the sample
    ## weight turns the ensemble into the flux-weighted one; the detector then sees a
    ## crossing flux, and dividing by |v_x,det| converts that flux back into f.
    vxi = [s.u[1][4] for s in sols.u]
    vs, ws_init = get_particle_crossings(sols, detector, vxi)

    v_edges = -1000:dv_km:1000
    centers = bin_centers(v_edges)
    ## Each macro-particle stands for `n0 / N` of the source density spread over one bin
    ## volume `dv³`, hence `S = n0 [km⁻³] / (N · dv_km³)` for `f` in [s³/km⁶].
    S = (n0 * 1.0e9) / (length(sols.u) * dv_km^3)

    f_3d = bin_velocity_space(vs, fill(S, length(vs)), v_edges; vx_source = ws_init)
    f_xy, f_xz, f_yz = project_vdf(f_3d, dv_km)

    return ((centers, centers, f_xy), (centers, centers, f_xz), (centers, centers, f_yz))
end

# ## Plotting helpers
#
# A single log-scale colour map is used everywhere.  `compute_common_cr` finds a global
# `(fmin, fmax)` across a set of 2-D distributions so that the comparison plot uses one
# shared colour range, making physical differences visible rather than hidden by per-panel
# scaling.

"""
Display a 2-D distribution given as an `(x, y, A)` tuple on a log10 colour scale.
Zeros / NaNs are clamped to `fmin` so empty cells use the bottom of the colour bar.
"""
function logheatmap!(ax, h::Tuple; colormap = :turbo, cr = (1.0e-8, 1.0))
    x, y, A = h
    fmin, fmax = cr
    B = [isfinite(v) && v ≥ fmin ? min(v, fmax) : fmin for v in A]
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
        m = maximum(h[3])
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
    titles = ["Monte Carlo", "Forward Liouville", "Backward Liouville"]
    xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
    ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]
    hists = (h1, h2, h3)
    all_h = (h1..., h2..., h3...)
    cr = compute_common_cr(all_h)
    fig = Figure(size = (1300, 1100), fontsize = 22)
    gl = fig[1, 1] = GridLayout()
    Label(
        gl[1, 2:4], "Downstream velocity distributions (x = $(x_downstream * 1.0e-3) km)";
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
                Colorbar(
                    gl[2:4, 5], hm;
                    label = L"\log_{10}([\mathrm{s}^2/\mathrm{km}^5])",
                    labelsize = 18, ticklabelsize = 14
                )
            end
        end
        Label(gl[row, 1], titles[r]; fontsize = 24, rotation = π / 2, tellheight = false)
    end
    return fig
end

hists_up = reconstruct_mc_projections(sols, detector_up, n_up, 20.0)
hists_down = reconstruct_mc_projections(sols, detector_down, n_up, 20.0)

fig_mc = plot_shock_vdf(hists_up, hists_down, x_upstream, x_downstream)
fig_mc = DisplayAs.PNG(fig_mc) #hide

# Each crossing contributes `S · |v_x,src| / |v_x,det|`, so the histogram estimates the
# phase-space density ``f`` at the detector: the `N` and `dv` factors are absorbed in `S`,
# and the ``|v_x|`` ratio converts the crossing flux into a density.
#
# ## Method 2: Forward Liouville Tracking
#
# Forward Liouville tracking starts from a sphere of initial conditions in velocity space at
# the source and traces forward to the detector.  By Liouville's theorem
# ``f_{\det}(\mathbf{v}_{\det}) = f_{\src}(\mathbf{v}_{\src})``; the source value
# ``n_0\,\mathrm{pdf}(\mathrm{vdf}, \mathbf{v}_{\src})`` is carried unchanged along each
# trajectory.  Because the sphere is sampled uniformly, every sample stands for a known
# source velocity volume ``V_{\sph}/N``, which the trajectory maps onto a detector volume
# ``(V_{\sph}/N)\,|v_{x,\src}|/|v_{x,\det}|``.  Depositing ``f\,\Delta V`` into the detector
# bin and dividing by the bin volume gives the same bin-averaged ``f`` as Method 3, without
# the sampling noise of Method 1.

function reconstruct_liouville_projections(
        sols, detector, vdf, n0;
        dv_km = 20.0, vsphere = (4 / 3) * π * (vradius_m2 * 1.0e-3)^3
    )
    ## Source f for every trajectory [s³/m⁶]
    ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols.u]
    vxi = [s.u[1][4] for s in sols.u]
    ## Detector crossings carry the source f (Liouville) together with the launch |v_x|
    vs, ws = get_particle_crossings(sols, detector, ws0)
    _, ws_vxi = get_particle_crossings(sols, detector, vxi)

    v_edges = -1000:dv_km:1000
    centers = bin_centers(v_edges)
    ## Each sample stands for a source velocity-space volume `vsphere / N`, which
    ## `bin_velocity_space` maps onto a detector volume `vsphere / N · |v_x,src|/|v_x,det|`
    ## and spreads over the bin volume, so the bin density is Σ f·ΔV / dv³. This summed
    ## estimator carries the `∝ 1/√n` counting noise of the samples per bin; the
    ## `average = true` ratio estimator divides it out instead, at the cost of giving every
    ## visited bin the local `f` regardless of how many times a trajectory crossed it.
    ## With N = 10⁵ samples over the sphere the summed estimator is populated well enough
    ## here, so we keep the more faithful rendering of the mapped velocity volume.
    scale = vsphere * 1.0e18 / (length(sols.u) * dv_km^3) # [s³/m⁶] → [s³/km⁶]

    f_3d = bin_velocity_space(vs, ws .* scale, v_edges; vx_source = ws_vxi)
    f_xy, f_xz, f_yz = project_vdf(f_3d, dv_km)

    return ((centers, centers, f_xy), (centers, centers, f_xz), (centers, centers, f_yz))
end

nparticles_m2 = 100000
const vradius_m2 = 3 * vth_ion # velocity-space radius, [m/s]

## Uniform sampling in a 3D sphere
function prob_func_m2(prob, ctx)
    v = sample_velocity_ball(ctx.rng, vradius_m2; center = SA[V_sw, 0.0, 0.0])
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

prob_m2 = TraceProblem(
    SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param; prob_func = prob_func_m2
)
t_liou = @elapsed sols_m2 = TP.solve(
    prob_m2, Boris(), EnsembleThreads(); dt, savestepinterval = 1,
    trajectories = nparticles_m2, seed
);

hists_up_m2 = reconstruct_liouville_projections(sols_m2, detector_up, vdf, n_up)
hists_down_m2 = reconstruct_liouville_projections(sols_m2, detector_down, vdf, n_up)

fig_forward = plot_shock_vdf(hists_up_m2, hists_down_m2, x_upstream, x_downstream)
fig_forward = DisplayAs.PNG(fig_forward) #hide

# The sphere radius `3 vth_ion` covers the bulk of the Maxwellian; `nparticles_m2 = 10⁵`
# gives usable statistics in the populated bins. Note that the sharp circular boundary in
# the reconstructed phase-space plots is an artifact of the finite sampling sphere
# (`r ≤ 3 vth_ion`) at the source. The result is a direct Monte-Carlo estimate of
# ``f_{\det}`` on the same grid and in the same units as Methods 1 and 3, so all three
# should agree up to sampling noise.
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

const source_plane = Meshes.Plane(Meshes.Point(x_source...), Meshes.Vec(1.0, 0.0, 0.0))

## One backward-tracing pass over a uniform velocity grid; returns the 3-D phase-space
## density at the detector in [s³/km⁶].
function run_backward_pass(vx_grid, vy_grid, vz_grid, detector_x, dt, param)
    ## 20 s backward span: setting post_source_margin = 100 km safely captures the
    ## post-crossing point.
    ## Trajectories moving deeper downstream (u[1] < detector_x - 600 km) cannot return
    ## through the shock ramp and are terminated early, yielding an order-of-magnitude
    ## speedup.
    ## Note on savestepinterval: With dt ≈ τ_g / 20, savestepinterval = 10 spaces saved
    ## points by 180° of gyro-phase (half an orbit); linear interpolation on a semicircle
    ## collapses the perpendicular velocity toward zero and distorts f ∝ exp(-v²/(2vth²)).
    ## We keep savestepinterval = 1 for accurate boundary interpolation.
    post_source_margin = 100.0e3
    prob = vdf_grid_problem(
        vx_grid, vy_grid, vz_grid, SA[detector_x, 0.0, 0.0], param, (0.0, -20.0)
    )

    sols = TP.solve(
        prob, Boris(), EnsembleThreads(); dt = -dt,
        trajectories = length(vx_grid) * length(vy_grid) * length(vz_grid),
        savestepinterval = 1,
        isoutside = (u, p, t) -> u[1] > x_source[1] + post_source_margin ||
            u[1] < detector_x - 600.0e3
    )

    return vdf_backward(
        sols, source_plane, f_src,
        (length(vx_grid), length(vy_grid), length(vz_grid))
    )
end

function reconstruct_backward_projections(
        detector_x, dt, param;
        v_range = 1000.0e3, vy_range = 1000.0e3, vz_range = 1000.0e3, dv_km = 20.0,
        adaptive = true, dv_coarse_km = 60.0, margin_km = 150.0
    )
    dv = dv_km * 1.0e3
    ## Grid points are bin *centres*, matching the histograms of Methods 1 & 2 whose bin
    ## edges run from -v_range to +v_range in steps of dv. Snapping the refined window to
    ## those centres keeps all three methods on the same velocity grid.
    v0x = -v_range + dv / 2
    v0y = -vy_range + dv / 2
    v0z = -vz_range + dv / 2

    if adaptive
        ## Pass 1 (coarse): locate the populated region cheaply.
        vx_c = range(-v_range, v_range, step = dv_coarse_km * 1.0e3)
        vy_c = range(-vy_range, vy_range, step = dv_coarse_km * 1.0e3)
        vz_c = range(-vz_range, vz_range, step = dv_coarse_km * 1.0e3)
    else
        vx_c = range(v0x, -v0x, step = dv)
        vy_c = range(v0y, -v0y, step = dv)
        vz_c = range(v0z, -v0z, step = dv)
    end

    t_solve = @elapsed begin
        f_coarse = run_backward_pass(vx_c, vy_c, vz_c, detector_x, dt, param)
        if adaptive
            vx_grid, vy_grid, vz_grid = refine_vdf_window(
                f_coarse, vx_c, vy_c, vz_c, (v0x, v0y, v0z), dv;
                margin = margin_km * 1.0e3, relthresh = 1.0e-6,
                bounds = ((v0x, -v0x), (v0y, -v0y), (v0z, -v0z))
            )
            f_3d_km = run_backward_pass(vx_grid, vy_grid, vz_grid, detector_x, dt, param)
        else
            vx_grid, vy_grid, vz_grid = vx_c, vy_c, vz_c
            f_3d_km = f_coarse
        end
    end
    nparticles_bw = length(vx_grid) * length(vy_grid) * length(vz_grid)

    ## Integral over third dimension [km/s]: f_int in [s²/km⁵]
    f_xy, f_xz, f_yz = project_vdf(f_3d_km, dv_km)

    ## Embed into full grid for seamless display matching Methods 1 & 2
    full_centers = collect(range(v0x, -v0x; step = dv) .* 1.0e-3)
    function embed_2d(g1, g2, M)
        full_M = zeros(length(full_centers), length(full_centers))
        i1 = round(Int, (g1[1] - full_centers[1]) / dv_km) + 1
        i2 = round(Int, (g2[1] - full_centers[1]) / dv_km) + 1
        i1_end = min(length(full_centers), i1 + length(g1) - 1)
        i2_end = min(length(full_centers), i2 + length(g2) - 1)
        len1 = i1_end - i1 + 1
        len2 = i2_end - i2 + 1
        full_M[i1:i1_end, i2:i2_end] .= @view M[1:len1, 1:len2]
        return (full_centers, full_centers, full_M)
    end

    return (
            embed_2d(vx_grid .* 1.0e-3, vy_grid .* 1.0e-3, f_xy),
            embed_2d(vx_grid .* 1.0e-3, vz_grid .* 1.0e-3, f_xz),
            embed_2d(vy_grid .* 1.0e-3, vz_grid .* 1.0e-3, f_yz),
        ), t_solve, nparticles_bw
end

res_up_bw, t_bw_up, n_bw_up = reconstruct_backward_projections(x_upstream, dt, param)
res_down_bw, t_bw_down, n_bw_down = reconstruct_backward_projections(x_downstream, dt, param)
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

# ## Moment check: density and momentum
#
# Integrating a 2-D projection over velocity returns the lowest moments of the reconstructed
# distribution: the density ``n = \int f\,\mathrm{d}v_i\mathrm{d}v_j`` and the momentum
# ``n\,V_i = \int v_i f\,\mathrm{d}v_i\mathrm{d}v_j``. Two different statements can be read
# off these numbers.
#
# 1. **The three methods must agree with each other.** They all reconstruct the same ``f`` in
#    the same units, so their moments have to match; any mismatch is an error in the
#    estimator rather than in the physics. Because the moments integrate the whole VDF, they
#    are sensitive to normalisation errors and to repeated counting of the same trajectory,
#    neither of which is obvious on a log colour scale.
# 2. **The upstream row is an absolute calibration.** The upstream detector sits in the
#    uniform, undisturbed solar wind, where the answer is known a priori:
#    ``n = n_{\up}`` and ``n\,V_x = n_{\up} V_{\sw}``. Matching those values is a genuine
#    validation, not merely a consistency check.
#
# What these numbers do *not* test is the downstream state. The fields here are prescribed
# rather than solved self-consistently, and the particles are launched as a single pulse, so
# mass flux need not be conserved across the ramp: part of the beam is turned back inside the
# ramp before it can reach either detector. The upstream to downstream deficit visible below
# is a property of the model, not a reconstruction error.

using Markdown, Printf #hide
io_m = IOBuffer() #hide
println(io_m, "| Method | n up [10⁶ m⁻³] | n·Vx up [10⁹ m⁻³·km/s] | n down [10⁶ m⁻³] | n·Vx down [10⁹ m⁻³·km/s] |") #hide
println(io_m, "| :--- | :---: | :---: | :---: | :---: |") #hide
for (i, name) in enumerate(("Monte Carlo", "Forward Liouville", "Backward Liouville")) #hide
    mu = velocity_moments(((hists_up, hists_up_m2, res_up_bw)[i])[1]; dv = 20.0) #hide
    md = velocity_moments(((hists_down, hists_down_m2, res_down_bw)[i])[1]; dv = 20.0) #hide
    @printf( #hide
        io_m, "| **%s** | %.2f | %.2f | %.2f | %.2f |\n", #hide
        name, mu.n * 1.0e-6, mu.nV * 1.0e-9, md.n * 1.0e-6, md.nV * 1.0e-9 #hide
    ) #hide
end #hide
Markdown.parse(String(take!(io_m))) #hide

# ## Bin-wise deviation
#
# The moments collapse each VDF into four numbers, so a wrong shape with the right
# normalisation and the right shape with the wrong normalisation can look identical to them.
# For a pointwise comparison we use the same relative L2 norm as the steady-state demo,
# ``\|f_{\rec} - f_{\ref}\| / \|f_{\ref}\|``, evaluated over the populated cells of the bin
# grid that all three methods now share.
#
# Upstream the reference can be exact: that plane still sees the undisturbed solar wind, and
# since ``\mathbf{E} = -\mathbf{V}_{\sw}\times\mathbf{B}`` there the drifting Maxwellian is an
# exact steady solution, so `vdf` itself is the reference. Downstream no analytic reference
# exists, and the noise-free backward solution takes that role instead.

## Analytic references come out as bare matrices, while the reconstructions are
## `(vi, vj, f)` tuples; this normalizes both to the same matrix.
matrix_of(h::Tuple) = h[3]
matrix_of(M::AbstractMatrix) = M

const bin_edges = -1000.0:20.0:1000.0 # km/s, velocity bin edges of all three methods
const v_centers = bin_centers(bin_edges)
## Midpoint rule on the bin centres, i.e. the same quadrature the reconstructions use when
## they collapse the third velocity axis.
const v_int = range(-990.0, 990.0; step = 20.0) # km/s, integration axis

ana_up = (
    analytic_projection(f_src, 3, v_centers, v_centers, v_int, step(v_int)),
    analytic_projection(f_src, 2, v_centers, v_centers, v_int, step(v_int)),
    analytic_projection(f_src, 1, v_centers, v_centers, v_int, step(v_int)),
)

io_d = IOBuffer() #hide
println(io_d, "| Comparison | Reference | Vx–Vy | Vx–Vz | Vy–Vz |") #hide
println(io_d, "| :--- | :--- | :---: | :---: | :---: |") #hide
for (name, rec, ref, refname) in ( #hide
        ("Monte Carlo (up)", hists_up, ana_up, "analytic solar wind"), #hide
        ("Forward Liouville (up)", hists_up_m2, ana_up, "analytic solar wind"), #hide
        ("Backward Liouville (up)", res_up_bw, ana_up, "analytic solar wind"), #hide
        ("Monte Carlo (down)", hists_down, res_down_bw, "Backward Liouville"), #hide
        ("Forward Liouville (down)", hists_down_m2, res_down_bw, "Backward Liouville"), #hide
    ) #hide
    @printf( #hide
        io_d, "| **%s** | %s | %.3f | %.3f | %.3f |\n", #hide
        name, refname, (relative_l2(matrix_of(rec[i]), matrix_of(ref[i])) for i in 1:3)... #hide
    ) #hide
end #hide
Markdown.parse(String(take!(io_d))) #hide

# ## Summary
#
# This example illustrates three complementary ways to reconstruct the phase-space density
# from particle simulations, all returning the same physical quantity ``f`` and sharing a
# common colour scale in the comparison plot.

t_per_mc = t_mc / nparticles * 1.0e6
t_per_liou = t_liou / nparticles_m2 * 1.0e6
t_per_bw = t_bw / n_bw * 1.0e6

io = IOBuffer() #hide
println(io, "| Method | Monte Carlo | Forward Liouville Tracking | Backward Liouville Tracing |") #hide
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
