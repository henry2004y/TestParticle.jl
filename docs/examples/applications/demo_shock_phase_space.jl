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
# In this method, simulated particles are treated as macro-particles.
# Instead of calculating the density by counting snapshots in time, we treat the source as a steady-state flux injection.
# To convert crossing events into physical phase-space density, we apply kinematic weighting that maps the instantaneous launch of macro-particles to a continuous stream in a steady state.

# ### Projections
# We bin the crossing events into 2D orthogonal velocity planes, integrating over the third dimension.

function bin_crossings_flux_projected(vs, ws, n0, trajectories; dv_km = 20.0)
    v_edges = -1000:dv_km:1000
    h_xy = Hist2D(; binedges = (v_edges, v_edges))
    h_xz = Hist2D(; binedges = (v_edges, v_edges))
    h_yz = Hist2D(; binedges = (v_edges, v_edges))

    ## Normalization S = n / (N * dv^2)
    dv2 = (dv_km * 1.0e3)^2
    S = n0 / (trajectories * dv2)

    for (v, vxi) in zip(vs, ws)
        ## Kinematic weight W = |vxi| / |vxd|
        ## This converts instantaneous launch counts to steady-state density.
        w = abs(vxi) / abs(v[1])

        vx_km, vy_km, vz_km = v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3
        push!(h_xy, vx_km, vy_km, w)
        push!(h_xz, vx_km, vz_km, w)
        push!(h_yz, vy_km, vz_km, w)
    end

    ## Units: [s^2/m^5] or similar. We scale by 1e9 for display.
    return h_xy * S, h_xz * S, h_yz * S
end

println("Calculating forward crossings (Flux Projections)...")
ws0 = [abs(s.u[1][4]) for s in sols] # initial |vx|
vs_up_flux, ws_up_flux = get_particle_crossings(sols, detector_up, ws0)
vs_down_flux, ws_down_flux = get_particle_crossings(sols, detector_down, ws0)

h_up_xy, h_up_xz, h_up_yz = bin_crossings_flux_projected(
    vs_up_flux, ws_up_flux, n0_p, trajectories
)
h_down_xy, h_down_xz, h_down_yz = bin_crossings_flux_projected(
    vs_down_flux, ws_down_flux, n0_p, trajectories
)

fig_flux = Figure(size = (1000, 600), fontsize = 18)
labels = [L"v_x-v_y", L"v_x-v_z", L"v_y-v_z"]
hists_up = [h_up_xy, h_up_xz, h_up_yz]
hists_down = [h_down_xy, h_down_xz, h_down_yz]

for i in 1:3
    ax_up = Axis(fig_flux[1, i], title = "Upstream $(labels[i])")
    hm_up = heatmap!(ax_up, hists_up[i], colormap = :turbo)
    if i == 3
        Colorbar(fig_flux[1, 4], hm_up, label = L"[\mathrm{s}^2/\mathrm{km}^5 \cdot 10^{-9}]")
    end

    ax_down = Axis(fig_flux[2, i], title = "Downstream $(labels[i])")
    hm_down = heatmap!(ax_down, hists_down[i], colormap = :turbo)
    if i == 3
        Colorbar(fig_flux[2, 4], hm_down, label = L"[\mathrm{s}^2/\mathrm{km}^5 \cdot 10^{-9}]")
    end
end
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
