# # Shock Phase Space
#
# This example demonstrates how to trace ions across a collisionless shock and analyze their phase space distribution, following the demo from IRF-matlab.
# We utilize Liouville's theorem (phase space density conservation) and backward/forward tracing to reconstruct the distribution function.
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

const Tio = 3.0  # eV
const vitho = sqrt(2 * TP.qᵢ * Tio / TP.mᵢ) # ~24.0 km/s

const Vsw = -400.0e3 # m/s
const nsw = 5.0e6  # m^-3
const P0_val = 0.08e-9 # Pa
const l_shock = 5.0e3 # m

## Shock Parameters
const n0_p = 3.0e6 # m^-3
const n1_p = 8.0e6 # m^-3

## Magnetic Field Parameters
const theta = 45.0
const Bmagnitude = 30.0e-9 # T

function get_B0_B1(theta, B_mag)
    theta_rad = deg2rad(theta)
    Bup_x = B_mag * cos(theta_rad)
    Bup_y = B_mag * sin(theta_rad)
    Bup_mag = B_mag

    Bdown_x = Bup_x
    Bdown_y = 3 * Bup_y
    Bdown_mag = sqrt(Bdown_x^2 + Bdown_y^2)

    B0 = 0.5 * (Bdown_mag - Bup_mag)
    B1 = 0.5 * (Bup_mag + Bdown_mag)
    return B0, B1
end

B0_calc, B1_calc = get_B0_B1(theta, Bmagnitude)
const B0_val = B0_calc # T
const B1_val = B1_calc # T
const Bx_val = 5.0e-9; # T

# ## Field Definitions
# We define custom analytical functions for the electric and magnetic fields across the shock.

"""
Magnetic Field
"""
function get_B_shock(r)
    x = r[1]
    by = -B0_val * tanh(x / l_shock) + B1_val
    bx = Bx_val
    bz = 0.0
    return SVector{3}(bx, by, bz)
end

"""
Electric Field based on generalized Ohm's law, including the Hall term and Electron Pressure Gradient.
"""
function get_E_shock(r)
    xnorm = r[1] / l_shock
    tanh_v = tanh(xnorm)
    sech_v = sech(xnorm)
    ## Ion density for Harris current sheet
    ni = -n0_p * tanh_v + n1_p

    ## Jz from Ampere's Law
    jz = -B0_val * sech_v^2 / (TP.μ₀ * l_shock)

    by = -B0_val * tanh_v + B1_val
    bx = Bx_val

    ## Ohm's Law and momentum equation terms
    eni = TP.qᵢ * ni
    ex = -jz * by / eni + P0_val * sech_v^2 / (eni * l_shock)
    ey = jz * bx / eni
    ez = -Vsw * (B1_val - B0_val)

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

# ## Analysis and Plotting
# We bin the particle trajectories into phase space histograms.

function bin_results(sols, n0, trajectories, dt_interp; dx_km = 5.0, dv_km = 10.0)
    println("Binning results...")
    xrange = [-300, 300] # km
    vrange = [-1000, 1000] # km/s

    x_edges = xrange[1]:dx_km:xrange[2]
    v_edges = vrange[1]:dv_km:vrange[2]

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
            x = state[1] / 1.0e3
            vx = state[4]
            vy = state[5]
            vz = state[6]

            vx_km = vx / 1.0e3
            vy_km = vy / 1.0e3
            vz_km = vz / 1.0e3

            ## Weight by local |vx| to get density from time-integrated counts
            w = abs(vx)
            push!(h_x_vx, x, vx_km, weight = w)
            push!(h_x_vy, x, vy_km, weight = w)
            push!(h_x_vz, x, vz_km, weight = w)
        end
    end
    println("Binning complete.")

    return h_x_vx * S, h_x_vy * S, h_x_vz * S
end

# Upstream density (consistent with VDF)
h_x_vx, h_x_vy, h_x_vz = bin_results(sols, n0_p, trajectories, dt_interp);

# ### Phase Space Plots
# The reconstructed phase space distributions.

fig = Figure(size = (800, 900), fontsize = 20)

ax1 = Axis(fig[1, 1], title = "Phase Space X-Vx", ylabel = "Vx [km/s]")
hm1 = heatmap!(ax1, h_x_vx, colormap = :turbo)
Colorbar(fig[1, 2], hm1, label = "f(x, vx) [s/m^4]")

ax2 = Axis(fig[2, 1], title = "Phase Space X-Vy", ylabel = "Vy [km/s]")
hm2 = heatmap!(ax2, h_x_vy, colormap = :turbo)
Colorbar(fig[2, 2], hm2, label = "f(x, vy) [s/m^4]")

ax3 = Axis(
    fig[3, 1], title = "Phase Space X-Vz", xlabel = "Position x [km]", ylabel = "Vz [km/s]"
)
hm3 = heatmap!(ax3, h_x_vz, colormap = :turbo)
Colorbar(fig[3, 2], hm3, label = "f(x, vz) [s/m^4]")
fig = DisplayAs.PNG(fig) #hide

# ### Backward Tracing and Reconstruction
# We now setup plane detectors at $x = 10^5$ and $x = -10^5$.
# We compare the forward reconstruction (from crossings) with the backward reconstruction at these locations.

## Define detectors
x_up, x_down = 1.0e5, -1.0e5
detector_up = Meshes.Plane(Meshes.Point(x_up, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))
detector_down = Meshes.Plane(Meshes.Point(x_down, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))

## Forward Reconstruction
println("Calculating forward crossings...")
vs_up, _ = get_particle_crossings(sols, detector_up)
vs_down, _ = get_particle_crossings(sols, detector_down)

## Backward Tracing Reconstruction
# Define a 2D velocity grid at the detector
vx_grid = range(-1000, 1000, length = 100) .* 1.0e3
vy_grid = range(-1000, 1000, length = 100) .* 1.0e3
vz_fixed = 0.0

function backward_trace(x_det, vx_grid, vy_grid, vz)
    f_det = fill(NaN, (length(vx_grid), length(vy_grid)))

    ## Backward tracing parameters
    tspan_bw = (0.0, -50.0)

    ## Use TraceProblem and Parallel Boris for backward tracing
    trajectories_bw = length(vx_grid) * length(vy_grid)
    u0_grid = [(vx, vy) for vx in vx_grid, vy in vy_grid]

    function prob_func_bw(prob, i, repeat)
        vx, vy = u0_grid[i]
        u0_bw = [x_det, 0.0, 0.0, vx, vy, vz]
        return remake(prob, u0 = u0_bw)
    end

    prob_bw = TraceProblem(
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan_bw, param;
        prob_func = prob_func_bw
    )

    println("Tracing $trajectories_bw points backward...")
    sols_bw = TP.solve(
        prob_bw, EnsembleThreads();
        dt = -dt_interp, trajectories = trajectories_bw, save_everystep = false
    )

    for i in 1:trajectories_bw
        sol = sols_bw[i]
        last_state = sol.u[end]
        if last_state[1] >= 300.0e3
            ## VDF at upstream: f(v_final)
            v_final = last_state[SA[4, 5, 6]]
            ## VDF(v) built-in probability assessment
            ## Standard 2D Maxwellian distribution f(vx, vy) integrated over vz
            ## Factor of sqrt(pi)*vth converts 3D PDF to 2D integrated PDF
            idx_x = (i - 1) % length(vx_grid) + 1
            idx_y = (i - 1) ÷ length(vx_grid) + 1
            f_det[idx_x, idx_y] =
                pdf(vdf, SA[v_final[1], v_final[2], 0.0]) * (sqrt(pi) * vitho)
        end
    end
    return f_det
end

println("Starting backward tracing at upstream detector...")
f_up_bw = backward_trace(x_up, vx_grid, vy_grid, vz_fixed)
println("Starting backward tracing at downstream detector...")
f_down_bw = backward_trace(x_down, vx_grid, vy_grid, vz_fixed)

## Visualization of Comparison
fig_comp = Figure(size = (1000, 800), fontsize = 20)

vx_grid_km = vx_grid ./ 1.0e3
vy_grid_km = vy_grid ./ 1.0e3

## Forward plots (Integrated f(vx, vy) at detector)
# To get density f from crossings, we weight by 1/|vx|
function bin_crossings(vs, n_upstream, trajectories; dv_km = 20.0)
    v_edges = -1000:dv_km:1000
    h = Hist2D(; binedges = (v_edges, v_edges))
    ## For a pulse/slab injection, the total counts of crossings N_cross(v)
    ## are directly proportional to the phase space density f(v).
    ## Normalization: S = n_upstream / (trajectories * dv^2)
    ## We multiply by 1e6 to convert s^2/m^2 to s^2/km^2.
    S = (n_upstream / trajectories) / (dv_km * 1.0e3)^2
    for v in vs
        push!(h, v[1] / 1.0e3, v[2] / 1.0e3)
    end
    return h * S * 1.0e6 # [s^2/km^2]
end

h_up_fw = bin_crossings(vs_up, n0_p, trajectories)
h_down_fw = bin_crossings(vs_down, n0_p, trajectories)

ax_up_fw = Axis(
    fig_comp[1, 1];
    title = "Forward f(vx, vy) at x = 100 km", xlabel = "vx [km/s]", ylabel = "vy [km/s]"
)
hm_up_fw = heatmap!(ax_up_fw, h_up_fw, colormap = :turbo)
Colorbar(fig_comp[1, 2], hm_up_fw, label = "f(vx, vy) [s^2/km^2]")

ax_down_fw = Axis(
    fig_comp[1, 3];
    title = "Forward f(vx, vy) at x = -100 km", xlabel = "vx [km/s]"
)
hm_down_fw = heatmap!(ax_down_fw, h_down_fw, colormap = :turbo)
Colorbar(fig_comp[1, 4], hm_down_fw, label = "f(vx, vy) [s^2/km^2]")

## Backward plots
ax_up_bw = Axis(
    fig_comp[2, 1];
    title = "Backward f(vx, vy) at x = 100 km", xlabel = "vx [km/s]", ylabel = "vy [km/s]"
)
hm_up_bw = heatmap!(
    ax_up_bw, vx_grid_km, vy_grid_km, f_up_bw .* 1.0e6;
    colormap = :turbo
)
Colorbar(fig_comp[2, 2], hm_up_bw, label = "f(vx, vy) [s^2/km^2]")

ax_down_bw = Axis(
    fig_comp[2, 3];
    title = "Backward f(vx, vy) at x = -100 km", xlabel = "vx [km/s]"
)
hm_down_bw = heatmap!(
    ax_down_bw, vx_grid_km, vy_grid_km, f_down_bw .* 1.0e6;
    colormap = :turbo
)
Colorbar(fig_comp[2, 4], hm_down_bw, label = "f(vx, vy) [s^2/km^2]")

fig_comp = DisplayAs.PNG(fig_comp) #hide
