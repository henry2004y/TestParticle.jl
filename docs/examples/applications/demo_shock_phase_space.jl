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

Random.seed!(42);

# ## Upstream Plasma Parameters

const T_ion = 20.0  # ion temperature [eV]
const vth_ion = sqrt(2 * TP.qᵢ * T_ion / TP.mᵢ) # ion thermal speed [m/s]
const V_sw = -400.0e3 # solar wind bulk speed [m/s]
const P_sw = 0.08e-9 # solar wind dynamic pressure [Pa]

# ## Shock Structure Parameters

const n_up = 3.0e6 # upstream number density [m⁻³]
const n_down = 8.0e6 # downstream number density [m⁻³]
const shock_width = 5.0e3 # shock ramp width [m]

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
const B_normal = 5.0e-9 # shock normal component of B [T]

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

const nparticles = 10000
const x_source = SA[300.0e3, 0.0, 0.0] # source plane location [m]
const tspan = (0.0, 20.0) # simulation time span [s]
const dt = get_gyroperiod(3 * B_mag) / 20 # time step [s]

param = prepare(get_E_shock, get_B_shock; species = Proton)

## Source velocity distribution (isotropic Maxwellian)
const p_thermal = n_up * TP.qᵢ * T_ion
const vdf = TP.Maxwellian(SA[V_sw, 0.0, 0.0], p_thermal, n_up; m = TP.mᵢ)

function prob_func_maxwellian(prob, i, repeat)
    v = rand(vdf)
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

u0_dummy = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func_maxwellian)

println("Starting simulation with $nparticles particles...")
@time sols = TP.solve(prob; dt, savestepinterval = 1, trajectories = nparticles);
println("Simulation complete.")

## Detector planes (upstream and downstream of the shock)
const x_upstream = 2.0e5  # [m]
const x_downstream = -2.0e5 # [m]
detector_up = Meshes.Plane(
    Meshes.Point(x_upstream, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0)
)
detector_down = Meshes.Plane(
    Meshes.Point(x_downstream, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0)
)


# ## Method 1: Flux Injection (The Macro-Particle Method)
# In this method, simulated particles are treated as macro-particles.
# Instead of calculating the density by counting snapshots in time, we treat the source as a steady-state flux injection.
# To convert crossing events into physical phase-space density, we apply kinematic weighting that maps the instantaneous launch of macro-particles to a continuous stream in a steady state.

# ### Projections
# We bin the crossing events into 2D orthogonal velocity planes, integrating over the third dimension.

function reconstruct_flux_projections(sols, detector, n0, dv_km)
    ## 1. Initial velocities at the source plane
    vxi = [s.u[1][4] for s in sols] # initial vx [m/s]
    ## 2. Detect crossings at the plane
    vs, ws_init = get_particle_crossings(sols, detector, vxi)

    v_edges = -1000:dv_km:1000
    h_3d = Hist3D(; binedges = (v_edges, v_edges, v_edges))

    ## Flux normalization factor S = n0_km3 / (N_total * dv_km^2)
    ## Conversion: n0 [m^-3] * 1e9 = n0 [km^-3]
    S = (n0 * 1.0e9) / (length(sols) * dv_km^2)

    for (v, vxi_val) in zip(vs, ws_init)
        ## Weight w = (v_xi / v_det) * S
        ## Resulting units: [s^2/km^5]
        w = abs(vxi_val) / abs(v[1]) * S
        push!(h_3d, v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3, w)
    end

    return project(h_3d, :z), project(h_3d, :y), project(h_3d, :x)
end

println("Calculating flux projections labels...")
h_up_xy, h_up_xz, h_up_yz = reconstruct_flux_projections(sols, detector_up, n_up, 20.0)
h_down_xy, h_down_xz, h_down_yz = reconstruct_flux_projections(sols, detector_down, n_up, 20.0)

fig_flux = Figure(size = (1200, 600), fontsize = 20)
hists_up = [h_up_xy, h_up_xz, h_up_yz]
hists_down = [h_down_xy, h_down_xz, h_down_yz]
xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]

for i in 1:3
    ax_up = Axis(
        fig_flux[1, i], title = "Upstream x = $(x_upstream * 1.0e-3) km",
        xlabel = xlabels[i], ylabel = ylabels[i]
    )
    hm_up = heatmap!(ax_up, hists_up[i], colormap = :turbo)
    if i == 3
        Colorbar(
            fig_flux[1, 4], hm_up;
            label = L"[\mathrm{s}^2/\mathrm{km}^5]"
        )
    end

    ax_down = Axis(
        fig_flux[2, i], title = "Downstream x = $(x_downstream * 1.0e-3) km",
        xlabel = xlabels[i], ylabel = ylabels[i]
    )
    hm_down = heatmap!(ax_down, hists_down[i], colormap = :turbo)
    if i == 3
        Colorbar(
            fig_flux[2, 4], hm_down;
            label = L"[\mathrm{s}^2/\mathrm{km}^5]"
        )
    end
end
fig_flux = DisplayAs.PNG(fig_flux) #hide

function reconstruct_liouville_projections(sols, detector, vdf, n0, Vsphere; dv_km = 20.0)
    ## 1. Initial weights from source PDF
    ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols]
    ## 2. Crossings
    vs, ws = get_particle_crossings(sols, detector, ws0)

    v_edges = -1000:dv_km:1000
    h_3d = Hist3D(; binedges = (v_edges, v_edges, v_edges))

    ## Normalization in km-based units: Vsphere [km^3], dv_km [km/s]
    ## Conversion: 1 m^3 = 1e-9 km^3
    Vsphere_km = Vsphere * 1.0e-9
    S_L = Vsphere_km / (length(sols) * dv_km^2)

    for (v, w) in zip(vs, ws)
        ## w is in [s^3/m^6]. Convert to [s^3/km^6] by multiplying 1e18.
        push!(h_3d, v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3, w * 1.0e18 * S_L)
    end

    return project(h_3d, :z), project(h_3d, :y), project(h_3d, :x)
end

println("Starting Liouville pathfinders...")
trajectories_m2 = 10000
const vradius_m2 = 3 * vth_ion # velocity space radius, [m/s]
const Vsphere_m2 = (4 / 3) * π * vradius_m2^3 # velocity space volume

## Uniform sampling in a 3D sphere
function prob_func_m2(prob, i, repeat)
    r = vradius_m2 * rand()^(1 / 3)
    ϕ = 2π * rand()
    θ = acos(2 * rand() - 1)

    sinθ, cosθ = sincos(θ)
    cosϕ, sinϕ = sincos(ϕ)
    v = SA[V_sw + r * sinθ * cosϕ, r * sinθ * sinϕ, r * cosθ]
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

prob_m2 = TraceProblem(SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param; prob_func = prob_func_m2)
@time sols_m2 = TP.solve(prob_m2; dt, savestepinterval = 1, trajectories = trajectories_m2);

println("Calculating forward crossings with conserved weights (Liouville Method)...")
h_up_xy_m2, h_up_xz_m2, h_up_yz_m2 = reconstruct_liouville_projections(
    sols_m2, detector_up, vdf, n_up, Vsphere_m2
)
h_down_xy_m2, h_down_xz_m2, h_down_yz_m2 = reconstruct_liouville_projections(
    sols_m2, detector_down, vdf, n_up, Vsphere_m2
)

fig_forward = Figure(size = (1200, 600), fontsize = 20)
hists_up_m2 = [h_up_xy_m2, h_up_xz_m2, h_up_yz_m2]
hists_down_m2 = [h_down_xy_m2, h_down_xz_m2, h_down_yz_m2]
xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]

for i in 1:3
    ax_up = Axis(
        fig_forward[1, i], title = "Upstream x = $(x_upstream * 1.0e-3) km",
        xlabel = xlabels[i], ylabel = ylabels[i]
    )
    ## No scaling for consistency with raw s^2/km^5 units
    hm_up = heatmap!(ax_up, hists_up_m2[i], colormap = :turbo)
    if i == 3
        Colorbar(fig_forward[1, 4], hm_up, label = L"[\mathrm{s}^2/\mathrm{km}^5]")
    end

    ax_down = Axis(
        fig_forward[2, i], title = "Downstream x = $(x_downstream * 1.0e-3) km",
        xlabel = xlabels[i], ylabel = ylabels[i]
    )
    hm_down = heatmap!(ax_down, hists_down_m2[i], colormap = :turbo)
    if i == 3
        Colorbar(fig_forward[2, 4], hm_down, label = L"[\mathrm{s}^2/\mathrm{km}^5]")
    end
end
fig_forward = DisplayAs.PNG(fig_forward) #hide

# ## Method 3: Backward Tracing
# In backward tracing, we start from a grid in velocity space at the detector and trace backward.
# The phase space density at the detector is simply the source density evaluated at the traced initial state.
# We then integrate the resulting 3D grid of values to provide a comparison for the other methods.

function reconstruct_backward_projections(
        detector_x, vdf, n0, dt, param;
        v_range = 1000.0e3, dv_km = 40.0
    )
    vx_grid = range(-v_range, v_range, step = dv_km * 1.0e3)
    vy_grid = range(-v_range, v_range, step = dv_km * 1.0e3)
    vz_grid = range(-500.0e3, 500.0e3, step = dv_km * 1.0e3)
    dvz_km = step(vz_grid) * 1.0e-3 # km/s

    nx, ny, nz = length(vx_grid), length(vy_grid), length(vz_grid)
    trajectories_bw = nx * ny * nz
    println("Tracing $trajectories_bw points backward for x = $(detector_x * 1.0e-3) km...")

    ## Initial conditions at detector
    function prob_func_bw(prob, i, repeat)
        iz = (i - 1) % nz + 1
        iy = ((i - 1) ÷ nz) % ny + 1
        ix = ((i - 1) ÷ (nz * ny)) % nx + 1

        u0_bw = SA[detector_x, 0.0, 0.0, vx_grid[ix], vy_grid[iy], vz_grid[iz]]
        return remake(prob, u0 = u0_bw)
    end

    tspan_bw = (0.0, -20.0)
    is_at_source(u, p, t) = u[1] >= x_source[1]
    prob_bw = TraceProblem(
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan_bw, param;
        prob_func = prob_func_bw
    )

    sols_bw = TP.solve(
        prob_bw, EnsembleThreads(); dt = -dt, trajectories = trajectories_bw,
        save_everystep = false, isoutofdomain = is_at_source
    )

    ## Evaluate PDF at source for each traced state
    f_3d_km = zeros(nx, ny, nz)
    for i in 1:trajectories_bw
        sol = sols_bw[i]
        last_state = sol.u[end]
        if last_state[1] >= x_source[1]
            iz = (i - 1) % nz + 1
            iy = ((i - 1) ÷ nz) % ny + 1
            ix = ((i - 1) ÷ (nz * ny)) % nx + 1
            ## Evaluate f [s^3/m^6] and convert to [s^3/km^6]
            f_3d_km[ix, iy, iz] = n0 * pdf(vdf, last_state[SA[4, 5, 6]]) * 1.0e18
        end
    end

    ## Integral over third dimension [km/s]: f_int in [s^2/km^5]
    f_xy = dropdims(sum(f_3d_km, dims = 3), dims = 3) .* dvz_km
    f_xz = dropdims(sum(f_3d_km, dims = 2), dims = 2) .* (step(vy_grid) * 1.0e-3)
    f_yz = dropdims(sum(f_3d_km, dims = 1), dims = 1) .* (step(vx_grid) * 1.0e-3)

    return (vx_grid .* 1.0e-3, vy_grid .* 1.0e-3, f_xy),
        (vx_grid .* 1.0e-3, vz_grid .* 1.0e-3, f_xz),
        (vy_grid .* 1.0e-3, vz_grid .* 1.0e-3, f_yz)
end

println("Calculating backward tracing ground truth...")
res_up_bw = reconstruct_backward_projections(x_upstream, vdf, n_up, dt, param)
res_down_bw = reconstruct_backward_projections(x_downstream, vdf, n_up, dt, param)

fig_backward = Figure(size = (1200, 600), fontsize = 20)
xlabels = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
ylabels = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]

for i in 1:3
    ax_up = Axis(
        fig_backward[1, i], title = "Upstream x = $(x_upstream * 1.0e-3) km",
        xlabel = xlabels[i], ylabel = ylabels[i]
    )
    ## No scaling for consistency with raw s^2/km^5 units
    hm_up = heatmap!(
        ax_up, res_up_bw[i][1], res_up_bw[i][2], res_up_bw[i][3],
        colormap = :turbo
    )
    if i == 3
        Colorbar(fig_backward[1, 4], hm_up, label = L"[\mathrm{s}^2/\mathrm{km}^5]")
    end

    ax_down = Axis(
        fig_backward[2, i], title = "Downstream x = $(x_downstream * 1.0e-3) km",
        xlabel = xlabels[i], ylabel = ylabels[i]
    )
    hm_down = heatmap!(
        ax_down, res_down_bw[i][1], res_down_bw[i][2], res_down_bw[i][3],
        colormap = :turbo
    )
    if i == 3
        Colorbar(fig_backward[2, 4], hm_down, label = L"[\mathrm{s}^2/\mathrm{km}^5]")
    end
end
fig_backward = DisplayAs.PNG(fig_backward) #hide

# ## Summary
# This example illustrates three complementary ways to reconstruct the phase space density from particle simulations.
#
# | Method | Flux Injection | Forward Liouville Tracking | Backward Liouville Tracing |
# |:---|:---|:---|:---|
# | **Noise** | Statistical ($\propto 1/\sqrt{N}$) | Low (analytical weights) | None (grid-based) |
# | **Coverage** | Source-sampled | Source-sampled | Target-sampled |
# | **Cost** | High | Medium | Low |
# | **Tail resolution** | Poor without large $N$ | Limited by sphere radius | Uniform across grid |
# | **Post-processing** | Binning + weighting | Binning + projection | PDF evaluation only |
