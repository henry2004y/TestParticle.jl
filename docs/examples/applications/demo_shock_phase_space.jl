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
using Statistics
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
            push!(h_x_vx, x, vx_km, w)
            push!(h_x_vy, x, vy_km, w)
            push!(h_x_vz, x, vz_km, w)
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

# ### Backward Tracing and Reconstruction
# We now setup plane detectors at $x = 10^5$ and $x = -10^5$.
# We compare the forward reconstruction (from crossings) with the backward reconstruction at these locations.

## Define detectors
x_up, x_down = 1.0e5, -1.0e5 # [m]
detector_up = Meshes.Plane(Meshes.Point(x_up, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))
detector_down = Meshes.Plane(Meshes.Point(x_down, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))

## Forward Reconstruction
println("Calculating forward crossings...")
ws0 = [pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols]
vs_up, ws_up = get_particle_crossings(sols, detector_up, ws0)
vs_down, ws_down = get_particle_crossings(sols, detector_down, ws0)

## Calculate source flux Gamma for Framework 2
## Gamma = n * <|vx|> from the analytical distribution
vx_mean = mean(abs(v[1]) for v in (rand(vdf) for _ in 1:100000))
gamma_up = n0_p * vx_mean

## Backward Tracing Reconstruction
# Define a 2D velocity grid at the detector
vx_grid = range(-1000.0e3, 1000.0e3, length = 100)
vy_grid = range(-1000.0e3, 1000.0e3, length = 100)
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
        u0_bw = SA[x_det, 0.0, 0.0, vx, vy, vz]
        return remake(prob, u0 = u0_bw)
    end

    prob_bw = TraceProblem(
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan_bw, param;
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
            f_det[idx_x, idx_y] = pdf(vdf, v_final) * (sqrt(pi) * vitho)
        end
    end
    return f_det
end

println("Starting backward tracing at upstream detector...")
f_up_bw = backward_trace(x_up, vx_grid, vy_grid, vz_fixed)
println("Starting backward tracing at downstream detector...")
f_down_bw = backward_trace(x_down, vx_grid, vy_grid, vz_fixed)

## Visualization of Comparison
# Forward tracing plots (Integrated f(vx, vy) at detector)

## Framework 1: Phase Space Tracking (The Liouville Method)
## We take the arithmetic mean of conserved analytical weights from the source.
## To match the 2D integrated density, we multiply by sqrt(2*pi)*sigma.
function bin_crossings_liouville(vs, ws, vitho; dv_km = 20.0)
    v_edges = -1000:dv_km:1000
    h_sum = Hist2D(; binedges = (v_edges, v_edges))
    h_cnt = Hist2D(; binedges = (v_edges, v_edges))

    for (v, w) in zip(vs, ws)
        push!(h_sum, v[1] / 1.0e3, v[2] / 1.0e3, w)
        push!(h_cnt, v[1] / 1.0e3, v[2] / 1.0e3, 1.0)
    end
    ## Mean weight in each bin
    f_mean = h_sum.bincounts ./ h_cnt.bincounts
    ## Replace NaNs (no particles) with 0
    f_mean[isnan.(f_mean)] .= 0.0

    ## Convert 3D PDF to 2D integrated PDF: factor of sqrt(pi)*vitho
    return h_sum, f_mean .* (sqrt(pi) * vitho) * 1.0e6
end

## Framework 2: Flux Injection (The Macro-Particle Method)
## We treat particles as macro-particles and account for physical flux.
## f = counts / (dt * dA * dv^3 * |vx|)
function bin_crossings_flux(vs, trajectories, gamma; dv_km = 20.0)
    v_edges = -1000:dv_km:1000
    h = Hist2D(; binedges = (v_edges, v_edges))

    ## Weight each crossing by 1/|vx| to recover density from flux
    for v in vs
        vx = abs(v[1])
        push!(h, v[1] / 1.0e3, v[2] / 1.0e3, 1.0 / vx)
    end

    ## Normalization: S = Gamma / (trajectories * dv^2)
    ## We multiply by 1e6 to convert s^2/m^2 to s^2/km^2.
    S = (gamma / trajectories) / (dv_km * 1.0e3)^2
    return h * S * 1.0e6 # [s^2/km^2]
end

h_up_fw_m1_cnt, f_up_fw_m1 = bin_crossings_liouville(vs_up, ws_up, vitho)
h_down_fw_m1_cnt, f_down_fw_m1 = bin_crossings_liouville(vs_down, ws_down, vitho)

h_up_fw_m2 = bin_crossings_flux(vs_up, trajectories, gamma_up)
h_down_fw_m2 = bin_crossings_flux(vs_down, trajectories, gamma_up)


fig_comp = Figure(size = (1000, 1000), fontsize = 20)

vx_grid_km = vx_grid ./ 1.0e3
vy_grid_km = vy_grid ./ 1.0e3

## Row 1: Forward Method 1 (Liouville)
ax_up_m1 = Axis(
    fig_comp[1, 1];
    title = "Liouville Upstream (x = 100 km)", ylabel = "vy [km/s]"
)
hm_up_m1 = heatmap!(ax_up_m1, h_up_fw_m1_cnt.binedges[1], h_up_fw_m1_cnt.binedges[2], f_up_fw_m1, colormap = :turbo)
Colorbar(fig_comp[1, 2], hm_up_m1, label = L"[\mathrm{s}^2/\mathrm{km}^2]")

ax_down_m1 = Axis(
    fig_comp[1, 3];
    title = "Liouville Downstream (x = -100 km)"
)
hm_down_m1 = heatmap!(ax_down_m1, h_down_fw_m1_cnt.binedges[1], h_down_fw_m1_cnt.binedges[2], f_down_fw_m1, colormap = :turbo)
Colorbar(fig_comp[1, 4], hm_down_m1, label = L"[\mathrm{s}^2/\mathrm{km}^2]")

## Row 2: Forward Method 2 (Flux)
ax_up_m2 = Axis(
    fig_comp[2, 1];
    title = "Flux Upstream (x = 100 km)", ylabel = "vy [km/s]"
)
hm_up_m2 = heatmap!(ax_up_m2, h_up_fw_m2, colormap = :turbo)
Colorbar(fig_comp[2, 2], hm_up_m2, label = L"[\mathrm{s}^2/\mathrm{km}^2]")

ax_down_m2 = Axis(
    fig_comp[2, 3];
    title = "Flux Downstream (x = -100 km)"
)
hm_down_m2 = heatmap!(ax_down_m2, h_down_fw_m2, colormap = :turbo)
Colorbar(fig_comp[2, 4], hm_down_m2, label = L"[\mathrm{s}^2/\mathrm{km}^2]")

## Row 3: Backward plots
ax_up_bw = Axis(
    fig_comp[3, 1];
    title = "Backward Upstream (x = 100 km)", xlabel = "vx [km/s]", ylabel = "vy [km/s]"
)
hm_up_bw = heatmap!(
    ax_up_bw, vx_grid_km, vy_grid_km, f_up_bw .* 1.0e6;
    colormap = :turbo
)
Colorbar(fig_comp[3, 2], hm_up_bw, label = L"[\mathrm{s}^2/\mathrm{km}^2]")

ax_down_bw = Axis(
    fig_comp[3, 3];
    title = "Backward Downstream (x = -100 km)", xlabel = "vx [km/s]"
)
hm_down_bw = heatmap!(
    ax_down_bw, vx_grid_km, vy_grid_km, f_down_bw .* 1.0e6;
    colormap = :turbo
)
Colorbar(fig_comp[3, 4], hm_down_bw, label = L"[\mathrm{s}^2/\mathrm{km}^2]")

fig_comp = DisplayAs.PNG(fig_comp) #hide
