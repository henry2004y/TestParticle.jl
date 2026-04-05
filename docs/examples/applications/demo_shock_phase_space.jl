# # Shock Phase Space
#
# This example demonstrates how to trace ions across a collisionless shock and analyze their phase space distribution, following the demo from IRF-matlab.
# We utilize Liouville's theorem (phase space density conservation), backward/forward tracing, and flux injection to reconstruct the distribution function.
# In this specific example, we trace particles forward from a Maxwellian distribution upstream to seeing their evolution downstream.

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

# ## Physics Constants & Parameters
# We use analytical profiles to represent the shock transition.

const Tio = 3.0  # [eV]
const vitho = sqrt(2 * TP.qᵢ * Tio / TP.mᵢ) # ~24.0 [km/s]

const Vsw = -400.0e3 # [m/s]
const Psw = 0.08e-9 # [Pa]
const shock_width = 5.0e3 # [m]

## Shock Parameters
const n0_p = 3.0e6 # [m^-3]
const n1_p = 8.0e6 # [m^-3]

## Magnetic Field Parameters
const θ = 45.0 # [degree]
const Bmag = 30.0e-9 # [T]

function get_B0_B1(θ, B_mag)
    Bup_y, Bup_x = B_mag .* sincosd(θ)
    Bup_mag = B_mag

    Bdown_x = Bup_x
    Bdown_y = 3 * Bup_y
    Bdown_mag = sqrt(Bdown_x^2 + Bdown_y^2)

    B0 = 0.5 * (Bdown_mag - Bup_mag)
    B1 = 0.5 * (Bup_mag + Bdown_mag)
    return B0, B1
end

const B0, B1 = get_B0_B1(θ, Bmag) # [T]
const Bx = 5.0e-9; # [T]

# ## Field Definitions
# We define custom analytical functions for the electric and magnetic fields across the shock.

"""
Magnetic Field
"""
function get_B_shock(r)
    x = r[1]
    by = -B0 * tanh(x / shock_width) + B1
    bx = Bx
    bz = 0.0
    return SVector{3}(bx, by, bz)
end

"""
Electric Field based on generalized Ohm's law, including the Hall term and Electron Pressure Gradient.
"""
function get_E_shock(r)
    xnorm = r[1] / shock_width
    tanh_v = tanh(xnorm)
    sech_v = sech(xnorm)
    ## Ion density for Harris current sheet
    ni = -n0_p * tanh_v + n1_p

    ## Jz from Ampere's Law
    jz = -B0 * sech_v^2 / (TP.μ₀ * shock_width)

    by = -B0 * tanh_v + B1
    bx = Bx

    ## Ohm's Law and momentum equation terms
    eni = TP.qᵢ * ni
    ex = -jz * by / eni + Psw * sech_v^2 / (eni * shock_width)
    ey = jz * bx / eni
    ez = -Vsw * (B1 - B0)

    return SVector{3}(ex, ey, ez)
end;

# ## Simulation Setup

trajectories = 5000 # number of particles
tspan = (0.0, 50.0) # s
dt_interp = 1.0e-3 # s

## Prepare the Hamiltonian system
param = prepare(get_E_shock, get_B_shock; species = Proton)

## Generate random initial velocities (Maxwellian)
# Calculate pressures from temperatures in eV
p_orig = n0_p * TP.qᵢ * Tio
const vdf = TP.Maxwellian(SA[Vsw, 0.0, 0.0], p_orig, n0_p; m = TP.mᵢ)

## Initial Conditions (Upstream)
function prob_func(prob, i, repeat)
    v = rand(vdf)

    x0_start = 300.0e3
    y0_start = 0.0
    z0_start = 0.0
    u0 = SA[x0_start, y0_start, z0_start, v[1], v[2], v[3]]

    return remake(prob, u0 = u0)
end

u0_dummy = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
prob = TraceProblem(u0_dummy, tspan, param; prob_func)

println("Starting simulation with $trajectories particles...")
@time sols = TP.solve(prob; dt = dt_interp, savestepinterval = 1, trajectories);
println("Simulation complete.")

## Define detectors
x_up, x_down = 1.0e5, -1.0e5 # [m]
detector_up = Meshes.Plane(Meshes.Point(x_up, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))
detector_down = Meshes.Plane(Meshes.Point(x_down, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))


# ## Method 1: Flux Injection (The Macro-Particle Method)
# In this method, simulated particles are treated as macro-particles. Instead of calculating the density directly by counting, we treat the source as a steady flux injection.
# To convert crossing counts $C_{bin}$ into absolute phase space density, we account for the physical source flux $\Gamma$, the number of numerical particles, and the spatial crossing velocity $|v_x|$ of the bin. Faster particles spend less time in a given spatial volume, which geometrically reduces their volumetric density relative to their flux.

# ### Volume Analysis
# We bin the particle trajectories into phase space histograms.

function bin_results(sols, n0, trajectories, dt_interp; dx_km = 5.0, dv_km = 10.0)
    x_edges = -300:dx_km:300
    v_edges = -1000:dv_km:1000

    h_x_vx = Hist2D(; binedges = (x_edges, v_edges))
    h_x_vy = Hist2D(; binedges = (x_edges, v_edges))
    h_x_vz = Hist2D(; binedges = (x_edges, v_edges))

    ## Normalization factor for phase space density f(x, v_i)
    ## n0 = integral f d3v. We integrate over 2 dimensions in each 2D plot.
    ## f(x, v_x) = \int f(x, vx, vy, vz) dvy dvz.
    ## Weight per time step should be |vx| to recover f from time-integrated counts.
    ## Normalization: S = n0 * dt / (trajectories * dx * dv)
    S = n0 * dt_interp / (trajectories * dx_km * 1.0e3 * dv_km * 1.0e3)

    ## Binning loop
    for i in eachindex(sols)
        s = sols[i]

        for state in s.u
            x = state[1] * 1.0e-3
            vx = state[4]
            vy = state[5]
            vz = state[6]

            vx_km = vx * 1.0e-3
            vy_km = vy * 1.0e-3
            vz_km = vz * 1.0e-3

            ## Weight by local |vx| to get density from time-integrated counts
            w = abs(vx)
            push!(h_x_vx, x, vx_km, w)
            push!(h_x_vy, x, vy_km, w)
            push!(h_x_vz, x, vz_km, w)
        end
    end

    return h_x_vx * S, h_x_vy * S, h_x_vz * S
end

h_x_vx, h_x_vy, h_x_vz = bin_results(sols, n0_p, trajectories, dt_interp);

# ### Phase Space Plots
# The reconstructed phase space distributions.

fig = Figure(size = (800, 900), fontsize = 20)
ax1 = Axis(fig[1, 1], title = "Phase Space X-Vx", ylabel = "Vx [km/s]")
hm1 = heatmap!(ax1, h_x_vx, colormap = :turbo)
Colorbar(fig[1, 2], hm1, label = L"[\mathrm{s}/\mathrm{m}^4]")
ax2 = Axis(fig[2, 1], title = "Phase Space X-Vy", ylabel = "Vy [km/s]")
hm2 = heatmap!(ax2, h_x_vy, colormap = :turbo)
Colorbar(fig[2, 2], hm2, label = L"[\mathrm{s}/\mathrm{m}^4]")
ax3 = Axis(
    fig[3, 1], title = "Phase Space X-Vz", xlabel = "Position x [km]", ylabel = "Vz [km/s]"
)
hm3 = heatmap!(ax3, h_x_vz, colormap = :turbo)
Colorbar(fig[3, 2], hm3, label = L"[\mathrm{s}/\mathrm{m}^4]")
fig = DisplayAs.PNG(fig) #hide

# ### Flux-Based Detector Reconstruction

function bin_crossings_flux(
        vs, n0, trajectories;
        dv_km = 20.0, vz_threshold = 10.0
    )
    v_edges = -1000:dv_km:1000
    h = Hist2D(; binedges = (v_edges, v_edges))

    vz_limit = vz_threshold * 1.0e3
    for v in vs
        if abs(v[3]) < vz_limit
            push!(h, v[1] / 1.0e3, v[2] / 1.0e3)
        end
    end

    ## f3D = n0 / (N * dvx * dvy * dvz) * count
    dv3 = (dv_km * 1.0e3)^2 * (2 * vz_limit)
    S = n0 / (trajectories * dv3)
    return h * S * 1.0e9
end

println("Calculating forward crossings (Flux Method)...")
vs_up_flux, _ = get_particle_crossings(sols, detector_up)
vs_down_flux, _ = get_particle_crossings(sols, detector_down)

h_up_fw_m1 = bin_crossings_flux(vs_up_flux, n0_p, trajectories)
h_down_fw_m1 = bin_crossings_flux(
    vs_down_flux, n0_p, trajectories
)

fig_flux = Figure(size = (1000, 400), fontsize = 18)
ax_up = Axis(fig_flux[1, 1], title = L"Flux Upstream ($v_z=0$)", ylabel = "vy [km/s]")
hm_up = heatmap!(ax_up, h_up_fw_m1, colormap = :turbo)
Colorbar(fig_flux[1, 2], hm_up, label = L"[\mathrm{s}^3/\mathrm{km}^3 \cdot 10^{-9}]")

ax_down = Axis(fig_flux[1, 3], title = L"Flux Downstream ($v_z=0$)", ylabel = "vy [km/s]")
hm_down = heatmap!(ax_down, h_down_fw_m1, colormap = :turbo)
Colorbar(fig_flux[1, 4], hm_down, label = L"[\mathrm{s}^3/\mathrm{km}^3 \cdot 10^{-9}]")
fig_flux = DisplayAs.PNG(fig_flux) #hide


# ## Method 2: Phase Space Tracking (The Liouville Method)
# This method uses Liouville's theorem (phase space density conservation). Simulated particles are used as pathfinders to map the analytical density from the source to the detector. For every individual particle $i$ launched with initial velocity $\mathbf{v}_{0,i}$, its exact analytical phase space density $w_i = f(\mathbf{v}_{0,i})$ is recorded as a conserved weight.
# At the detector, we take the arithmetic mean of the weights of all particles that fall into a specific velocity bin to find the phase space density.

function bin_crossings_liouville(vs, ws; dv_km = 20.0, vz_threshold = 10.0)
    v_edges = -1000:dv_km:1000
    h_sum = Hist2D(; binedges = (v_edges, v_edges))
    h_cnt = Hist2D(; binedges = (v_edges, v_edges))

    vz_limit = vz_threshold * 1.0e3
    for (v, w) in zip(vs, ws)
        if abs(v[3]) < vz_limit
            push!(h_sum, v[1] / 1.0e3, v[2] / 1.0e3, w)
            push!(h_cnt, v[1] / 1.0e3, v[2] / 1.0e3, 1.0)
        end
    end
    ## Mean weight in each bin (represents f3D at vz=0)
    f_mean = h_sum.bincounts ./ h_cnt.bincounts
    ## Replace NaNs (no particles) with 0
    f_mean[isnan.(f_mean)] .= 0.0

    return h_sum, f_mean * 1.0e9 # [s^3/(km^3 * m^3)]
end

println("Calculating forward crossings with conserved weights (Liouville Method)...")
ws0 = [n0_p * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols]
vs_up, ws_up = get_particle_crossings(sols, detector_up, ws0)
vs_down, ws_down = get_particle_crossings(sols, detector_down, ws0)

h_up_fw_m2_cnt, f_up_fw_m2 = bin_crossings_liouville(vs_up, ws_up)
h_down_fw_m2_cnt, f_down_fw_m2 = bin_crossings_liouville(vs_down, ws_down)

fig_liouville = Figure(size = (1000, 400), fontsize = 20)
ax_up = Axis(fig_liouville[1, 1], title = L"Liouville Upstream ($v_z=0$)", ylabel = "vy [km/s]")
hm_up = heatmap!(ax_up, h_up_fw_m2_cnt.binedges[1], h_up_fw_m2_cnt.binedges[2], f_up_fw_m2, colormap = :turbo)
Colorbar(fig_liouville[1, 2], hm_up, label = L"[\mathrm{s}^3/\mathrm{km}^3 \cdot 10^{-9}]")

ax_down = Axis(fig_liouville[1, 3], title = L"Liouville Downstream ($v_z=0$)", ylabel = "vy [km/s]")
hm_down = heatmap!(ax_down, h_down_fw_m2_cnt.binedges[1], h_down_fw_m2_cnt.binedges[2], f_down_fw_m2, colormap = :turbo)
Colorbar(fig_liouville[1, 4], hm_down, label = L"[\mathrm{s}^3/\mathrm{km}^3 \cdot 10^{-9}]")
fig_liouville = DisplayAs.PNG(fig_liouville) #hide


# ## Method 3: Backward Tracing (The Ground Truth)
# In backward tracing, we start from a grid in velocity space at the detector location and trace particles backward in time to the source. If a particle reaches the source, the phase space density at the detector is the source density at the initial $(\mathbf{x}_0, \mathbf{v}_0)$.

# Define a 2D velocity grid at the detector
vx_grid = range(-1000.0e3, 1000.0e3, length = 100)
vy_grid = range(-1000.0e3, 1000.0e3, length = 100)
vz_fixed = 0.0

function backward_trace(x_det, vx_grid, vy_grid, vz)
    f_det = zeros((length(vx_grid), length(vy_grid)))

    ## Backward tracing parameters
    tspan_bw = (0.0, -50.0)

    ## Use TraceProblem and Parallel Boris for backward tracing
    trajectories_bw = length(vx_grid) * length(vy_grid)
    u0_grid = [(vx, vy) for vx in vx_grid, vy in vy_grid]

    function prob_func_bw(prob, i, repeat)
        vx, vy = u0_grid[i]
        u0_bw = SA[x_det, 0.0, 0.0, vx, vy, vz]
        return remake(prob, u0 = u0_bw)
    end

    is_at_source(u, p, t) = u[1] >= 300.0e3

    prob_bw = TraceProblem(
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan_bw, param;
        prob_func = prob_func_bw
    )

    println("Tracing $trajectories_bw points backward...")
    sols_bw = TP.solve(
        prob_bw, EnsembleThreads();
        dt = -dt_interp, trajectories = trajectories_bw, save_everystep = false,
        isoutofdomain = is_at_source
    )

    for i in 1:trajectories_bw
        sol = sols_bw[i]
        last_state = sol.u[end]
        if last_state[1] >= 300.0e3
            v_final = last_state[SA[4, 5, 6]]
            idx_x = (i - 1) % length(vx_grid) + 1
            idx_y = (i - 1) ÷ length(vx_grid) + 1
            f_det[idx_x, idx_y] = n0_p * pdf(vdf, v_final)
        end
    end
    return f_det
end

println("Starting backward tracing at upstream detector...")
f_up_bw = backward_trace(x_up, vx_grid, vy_grid, vz_fixed)
println("Starting backward tracing at downstream detector...")
f_down_bw = backward_trace(x_down, vx_grid, vy_grid, vz_fixed)

fig_backward = Figure(size = (1000, 400), fontsize = 18)
vx_grid_km = vx_grid ./ 1.0e3
vy_grid_km = vy_grid ./ 1.0e3

ax_up = Axis(fig_backward[1, 1], title = L"Backward Upstream ($v_z=0$)", xlabel = "vx [km/s]", ylabel = "vy [km/s]")
hm_up = heatmap!(ax_up, vx_grid_km, vy_grid_km, f_up_bw .* 1.0e9, colormap = :turbo)
Colorbar(fig_backward[1, 2], hm_up, label = L"[\mathrm{s}^3/\mathrm{km}^3 \cdot 10^{-9}]")

ax_down = Axis(fig_backward[1, 3], title = L"Backward Downstream ($v_z=0$)", xlabel = "vx [km/s]", ylabel = "vy [km/s]")
hm_down = heatmap!(ax_down, vx_grid_km, vy_grid_km, f_down_bw .* 1.0e9, colormap = :turbo)
Colorbar(fig_backward[1, 4], hm_down, label = L"[\mathrm{s}^3/\mathrm{km}^3 \cdot 10^{-9}]")
fig_backward = DisplayAs.PNG(fig_backward) #hide


# ## Summary
# This example illustrates three complementary ways to reconstruct phase space density from particle simulations.
# **Flux Injection** is simple and robust for macro-particle counting.
# **Liouville Tracking** preserves analytical density information and is less noisy in sparse regions.
# **Backward Tracing** provides the direct mapping from the detector to the source, serves as the ground truth validation.
