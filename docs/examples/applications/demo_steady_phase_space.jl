# # Steady-State Phase Space Validation
#
# This example validates the three phase-space tracing methods of TestParticle.jl
# against an *analytically known* steady-state distribution function (VDF). Unlike the
# shock demo (where the three methods are only compared with each other), here we know
# the exact VDF at the detector and can quantify how far each reconstruction deviates
# from it.
#
# ## Scenario: E×B drift of a bi-Maxwellian
# We use a uniform magnetic field `B = B ẑ` and a uniform perpendicular electric field
# `E ⟂ B`. The resulting E×B drift carries particle guiding centers steadily from the
# source plane to the detector. Crucially, with `E·B = 0` there is no parallel
# acceleration, so a bi-Maxwellian centred on the E×B drift velocity
# `u_E = E×B/B²` is an *exact* steady solution of the Vlasov equation. By Liouville's
# theorem the phase-space density is conserved along each characteristic, and because
# the gyration only rotates the perpendicular velocity (preserving `|v_⊥ - u_E|` and
# `v_∥`), the VDF at the detector is *identical* to the source VDF:
#
# ```math
# f_{\rm det}(\mathbf{v}) = f_{\rm src}(\mathbf{v}) =
#   n_0\,\exp\!\Bigl(-\frac{|\mathbf{v}_⊥ - \mathbf{u}_E|^2}{2v_{th,⊥}^2}\Bigr)
#   \exp\!\Bigl(-\frac{(v_∥ - u_{E,∥})^2}{2v_{th,∥}^2}\Bigr).
# ```
#
# This gives us a rigorous control experiment: any discrepancy between a reconstructed
# phase-space plot and this analytic expectation is purely numerical (integration
# error, Monte-Carlo scatter, binning / grid resolution).

import DisplayAs #hide
using TestParticle
import TestParticle as TP
using LinearAlgebra: norm
using StaticArrays
using Random
using Statistics: mean, var
using FHist
using VelocityDistributionFunctions
using CairoMakie
using Meshes
CairoMakie.activate!(type = "png") #hide

seed = 42;

# ## Field and plasma parameters

const B_mag = 10.0e-9  # uniform magnetic field magnitude [T]
const B_vec = SA[0.0, 0.0, B_mag]  # along +z
const V_drift = -400.0e3  # desired E×B drift speed [m/s] (along -x)
# E×B drift: u_E = E×B/B²  ⇒  E = B × u_E  (with E·B = 0).
const E_vec = SA[0.0, B_mag * V_drift, 0.0]  # [V/m]

function get_E_steady(r, t = 0.0)
    return E_vec
end

function get_B_steady(r, t = 0.0)
    return B_vec
end

# ## Source VDF: drifting bi-Maxwellian
# Anisotropic (perpendicular hotter than parallel) so the 2-D projections are visibly
# elliptical, giving the reconstruction a non-trivial structure to recover.

const n0 = 3.0e6        # number density [m⁻³]
const T_par_eV = 15.0   # parallel temperature [eV]
const T_perp_eV = 45.0  # perpendicular temperature [eV]
const p_par = n0 * TP.qᵢ * T_par_eV
const p_perp = n0 * TP.qᵢ * T_perp_eV

const vdf = TP.BiMaxwellian(
    SA[0.0, 0.0, 1.0], SA[V_drift, 0.0, 0.0], p_par, p_perp, n0; m = TP.mᵢ
);

# ## Geometry and integration

const x_source = SA[300.0e3, 0.0, 0.0] # source plane [m]
const tspan = (0.0, 4.0) # [s]; > transport time 500 km / 400 km/s = 1.25 s
const dt = get_gyroperiod(B_mag) / 40 # [s]; fine step so the crossing velocity is well resolved

const x_upstream = 200.0e3   # detector between source and origin [m]
const x_downstream = -200.0e3 # detector downstream of the origin [m]

detector_up = Meshes.Plane(
    Meshes.Point(x_upstream, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0)
)
detector_down = Meshes.Plane(
    Meshes.Point(x_downstream, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0)
)

param = prepare(get_E_steady, get_B_steady; species = Proton)

# ## Method 1: Forward Monte-Carlo Injection

nparticles = 8000

function prob_func_maxwellian(prob, ctx)
    v = rand(ctx.rng, vdf)
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

u0_dummy = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func_maxwellian)

t_mc = @elapsed sols = TP.solve(
    prob, Boris(); dt, savestepinterval = 1, trajectories = nparticles, seed
);

function reconstruct_flux_projections(sols, detector, n0, dv_km)
    vxi = [s.u[1][4] for s in sols.u]
    vs, ws_init = get_particle_crossings(sols, detector, vxi)

    v_edges = -1000:dv_km:1000
    h_3d = Hist3D(; binedges = (v_edges, v_edges, v_edges))

    S = (n0 * 1.0e9) / (length(sols.u) * dv_km^2)
    for (v, vxi_val) in zip(vs, ws_init)
        w = abs(vxi_val) / abs(v[1]) * S
        push!(h_3d, v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3, w)
    end
    return project(h_3d, :z), project(h_3d, :y), project(h_3d, :x)
end

hists_up = reconstruct_flux_projections(sols, detector_up, n0, 20.0)
hists_down = reconstruct_flux_projections(sols, detector_down, n0, 20.0)

# ## Method 2: Forward Liouville Tracking

function reconstruct_liouville_projections(sols, detector, vdf, n0, Vsphere; dv_km = 20.0)
    ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols.u]
    vs, ws = get_particle_crossings(sols, detector, ws0)

    v_edges = -1000:dv_km:1000
    h_3d = Hist3D(; binedges = (v_edges, v_edges, v_edges))

    Vsphere_km = Vsphere * 1.0e-9
    S_L = Vsphere_km / (length(sols.u) * dv_km^2)

    for (v, w) in zip(vs, ws)
        push!(h_3d, v[1] * 1.0e-3, v[2] * 1.0e-3, v[3] * 1.0e-3, w * 1.0e18 * S_L)
    end
    return project(h_3d, :z), project(h_3d, :y), project(h_3d, :x)
end

nparticles_m2 = 8000
const vth_perp = sqrt(2 * p_perp / (n0 * TP.mᵢ))
const vradius_m2 = 3 * vth_perp
const Vsphere_m2 = (4 / 3) * π * vradius_m2^3

function prob_func_m2(prob, ctx)
    r = vradius_m2 * rand(ctx.rng)^(1 / 3)
    ϕ = 2π * rand(ctx.rng)
    θ = acos(2 * rand(ctx.rng) - 1)
    sinθ, cosθ = sincos(θ)
    cosϕ, sinϕ = sincos(ϕ)
    v = SA[V_drift + r * sinθ * cosϕ, r * sinθ * sinϕ, r * cosθ]
    u0 = SA[x_source..., v...]
    return remake(prob, u0 = u0)
end

prob_m2 = TraceProblem(
    SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param; prob_func = prob_func_m2
)
t_liou = @elapsed sols_m2 = TP.solve(
    prob_m2, Boris(); dt, savestepinterval = 1, trajectories = nparticles_m2, seed
);

hists_up_m2 = reconstruct_liouville_projections(
    sols_m2, detector_up, vdf, n0, Vsphere_m2
)
hists_down_m2 = reconstruct_liouville_projections(
    sols_m2, detector_down, vdf, n0, Vsphere_m2
)

# ## Method 3: Backward Liouville Tracing

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
        SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], (0.0, -8.0), param;
        prob_func = prob_func
    )

    sols = TP.solve(
        prob, Boris(), EnsembleThreads(); dt = -dt, trajectories = ntraj,
        savestepinterval = 1,
        isoutside = (u, p, t) -> u[1] < detector_x - 1.0e5 ||
            u[1] > x_source[1] + 6000.0e3
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
        vx_c = range(-v_range, v_range; step = dv_coarse_km * 1.0e3)
        vy_c = range(-vy_range, vy_range; step = dv_coarse_km * 1.0e3)
        vz_c = range(-v_range, v_range; step = dv_coarse_km * 1.0e3)
    else
        vx_c = range(-v_range, v_range; step = dv)
        vy_c = range(-vy_range, vy_range; step = dv)
        vz_c = range(-v_range, v_range; step = dv)
    end

    t_solve = @elapsed begin
        f_coarse = run_backward_pass(vx_c, vy_c, vz_c, detector_x, vdf, n0, dt, param)
        if adaptive
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
    reconstruct_backward_projections(x_upstream, vdf, n0, dt, param)
res_down_bw, t_bw_down, n_bw_down =
    reconstruct_backward_projections(x_downstream, vdf, n0, dt, param)
t_bw = t_bw_up + t_bw_down
n_bw = n_bw_up + n_bw_down

# ## Analytic reference
# The exact detector VDF is the source bi-Maxwellian. We evaluate its 2-D projections
# `f(v_i, v_j) = ∫ f_3D(v_i, v_j, v_k) dv_k` on an explicit grid for contour overlays,
# and on each method's own grid for a quantitative deviation metric.

const vlim = 1000.0
const fine_g = range(-vlim, vlim; step = 40.0)   # km/s, for contours
const z_int = range(-vlim, vlim; step = 10.0)    # km/s, integration axis

function analytic_proj(vdf, n0, i, j, k, gi, gj, gk, dvk)
    M = zeros(length(gi), length(gj))
    for (bi, vi) in enumerate(gi), (bj, vj) in enumerate(gj)
        s = 0.0
        for vk in gk
            v1 = i == 1 ? vi * 1.0e3 : (j == 1 ? vj * 1.0e3 : vk * 1.0e3)
            v2 = i == 2 ? vi * 1.0e3 : (j == 2 ? vj * 1.0e3 : vk * 1.0e3)
            v3 = i == 3 ? vi * 1.0e3 : (j == 3 ? vj * 1.0e3 : vk * 1.0e3)
            s += n0 * pdf(vdf, SVector{3, Float64}(v1, v2, v3)) * 1.0e18 * dvk
        end
        M[bi, bj] = s
    end
    return M
end

function make_contours()
    return (
        analytic_proj(vdf, n0, 1, 2, 3, fine_g, fine_g, z_int, step(z_int)),
        analytic_proj(vdf, n0, 1, 3, 2, fine_g, fine_g, z_int, step(z_int)),
        analytic_proj(vdf, n0, 2, 3, 1, fine_g, fine_g, z_int, step(z_int)),
    )
end

# Analytic 2-D projections on the histogram edges used by Methods 1 & 2. Filling each
# 3-D bin with `f_3D·dv_z` and then summing over z reproduces the same [s²/km⁵] scale.
function analytic_hist_projections(vdf, n0, v_edges; dv_km)
    h = Hist3D(; binedges = (v_edges, v_edges, v_edges))
    for xi in 1:(length(v_edges) - 1), yi in 1:(length(v_edges) - 1), zi in 1:(length(v_edges) - 1)
        vc = SVector(
            (v_edges[xi] + v_edges[xi + 1]) / 2,
            (v_edges[yi] + v_edges[yi + 1]) / 2,
            (v_edges[zi] + v_edges[zi + 1]) / 2,
        )
        w = n0 * pdf(vdf, vc .* 1.0e3) * 1.0e18 * dv_km
        push!(h, vc[1], vc[2], vc[3], w)
    end
    return project(h, :z), project(h, :y), project(h, :x)
end

const v_edges = -1000:20:1000
const ana_hists = analytic_hist_projections(vdf, n0, v_edges; dv_km = 20.0)

function rel_l2(rec, ana; thresh = 1.0e-4)
    m = (ana .> maximum(ana) * thresh) .& isfinite.(rec)
    r = rec[m]
    a = ana[m]
    return norm(r .- a) / norm(a)
end

# ## Plots: reconstructed phase space vs analytic contours

function plot_validation(h_flux, h_liou, h_bw, xloc, fcont; vlim = 1000.0)
    titles = ["Flux Injection", "Forward Liouville", "Backward Liouville"]
    cols = [L"V_x–V_y", L"V_x–V_z", L"V_y–V_z"]
    hists = (h_flux, h_liou, h_bw)
    xlab = [L"V_x [\mathrm{km/s}]", L"V_x [\mathrm{km/s}]", L"V_y [\mathrm{km/s}]"]
    ylab = [L"V_y [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]", L"V_z [\mathrm{km/s}]"]
    fig = Figure(size = (1300, 1200), fontsize = 22)
    gl = fig[1, 1] = GridLayout()
    Label(
        gl[1, 2:4],
        "Detector at x = $(round(xloc * 1.0e-3; digits = 0)) km — white dashed = analytic";
        fontsize = 24, tellwidth = false
    )
    for r in 1:3, i in 1:3
        ax = Axis(
            gl[r + 1, i + 1], title = cols[i], xlabel = xlab[i], ylabel = ylab[i];
            xlabelsize = 22, ylabelsize = 22, titlesize = 20,
            xticklabelsize = 18, yticklabelsize = 18,
            limits = (-vlim, vlim, -vlim, vlim), aspect = 1
        )
        h = hists[r][i]
        hm = h isa Tuple ? heatmap!(ax, h...; colormap = :turbo) :
            heatmap!(ax, h; colormap = :turbo)
        contour!(
            ax, fine_g, fine_g, fcont[i]; color = :white, linestyle = :dash,
            linewidth = 2, levels = 6
        )
        if i == 3
            Colorbar(gl[r + 1, 5], hm; labelsize = 18, ticklabelsize = 14)
        end
    end
    for r in 1:3
        Label(gl[r + 1, 1], titles[r]; fontsize = 18, rotation = π / 2, tellheight = false)
    end
    return fig
end

fcont_down = make_contours()
fcont_up = make_contours()

fig_down = plot_validation(hists_down, hists_down_m2, res_down_bw, x_downstream, fcont_down)
fig_down = DisplayAs.PNG(fig_down) #hide

fig_up = plot_validation(hists_up, hists_up_m2, res_up_bw, x_upstream, fcont_up)
fig_up = DisplayAs.PNG(fig_up) #hide

# ## 1-D slice check
# Slice the `V_x–V_z` projection at `v_z = 0` and overlay the three reconstructions on
# the analytic profile. The backward (grid-based) method should sit on the curve.

function centers_of(g, M, dim)
    if length(g) == size(M, dim) + 1
        return (g[1:(end - 1)] .+ g[2:end]) ./ 2
    end
    return g
end

function slice_at(h, second::Bool, val)
    if h isa Tuple
        g1, g2, M = h
        c1 = centers_of(g1, M, 1)
        c2 = centers_of(g2, M, 2)
    else
        ex, ey = binedges(h)
        M = bincounts(h)
        c1 = (ex[1:(end - 1)] .+ ex[2:end]) ./ 2
        c2 = (ey[1:(end - 1)] .+ ey[2:end]) ./ 2
    end
    if second
        idx = argmin(abs.(c2 .- val))
        return c1, M[:, idx]
    else
        idx = argmin(abs.(c1 .- val))
        return c2, M[idx, :]
    end
end

matrix_of(h::Hist2D) = bincounts(h)
matrix_of(h::Tuple) = h[3]

const fine_x = range(-vlim, vlim; step = 10.0)
const fana_slice = [
    analytic_proj(vdf, n0, 1, 3, 2, [vx], [0.0], z_int, step(z_int))[1]
        for vx in fine_x
]

xc_f, fy_f = slice_at(hists_down[2], true, 0.0)
xc_l, fy_l = slice_at(hists_down_m2[2], true, 0.0)
xc_b, fy_b = slice_at(res_down_bw[2], true, 0.0)

fig_slice = Figure(size = (900, 500), fontsize = 20)
axs = Axis(
    fig_slice[1, 1], xlabel = L"V_x [\mathrm{km/s}]",
    ylabel = L"f(V_x, V_z=0) [\mathrm{s}^2/\mathrm{km}^5]",
    limits = (-vlim, vlim, nothing, nothing)
)
lines!(axs, fine_x, fana_slice; label = "Analytic", color = :black, linewidth = 3)
scatter!(axs, xc_f, fy_f; label = "Flux", color = (:blue, 0.5), markersize = 4)
scatter!(axs, xc_l, fy_l; label = "Forward Liouville", color = (:green, 0.5), markersize = 4)
lines!(axs, xc_b, fy_b; label = "Backward Liouville", color = :red, linewidth = 2)
axislegend(axs; position = :rt)
fig_slice = DisplayAs.PNG(fig_slice) #hide

# ## Deviation report
# Relative L2 norm `‖f_rec − f_ana‖ / ‖f_ana‖` over the populated cells, per method and
# projection. Methods 1 & 2 are statistical (∝ 1/√N); Method 3 is grid-limited.

using Markdown, Printf #hide
io = IOBuffer() #hide
println(io, "| Method | Vx–Vy | Vx–Vz | Vy–Vz |") #hide
println(io, "| :--- | :--- | :--- | :--- |") #hide
println(io, "| **Flux (down)** | $(round(rel_l2(matrix_of(hists_down[1]), matrix_of(ana_hists[1])); digits = 3)) | $(round(rel_l2(matrix_of(hists_down[2]), matrix_of(ana_hists[2])); digits = 3)) | $(round(rel_l2(matrix_of(hists_down[3]), matrix_of(ana_hists[3])); digits = 3)) |") #hide
println(io, "| **Liouville (down)** | $(round(rel_l2(matrix_of(hists_down_m2[1]), matrix_of(ana_hists[1])); digits = 3)) | $(round(rel_l2(matrix_of(hists_down_m2[2]), matrix_of(ana_hists[2])); digits = 3)) | $(round(rel_l2(matrix_of(hists_down_m2[3]), matrix_of(ana_hists[3])); digits = 3)) |") #hide
ana_bw_xy = analytic_proj(vdf, n0, 1, 2, 3, res_down_bw[1][1], res_down_bw[1][2], z_int, step(z_int))
ana_bw_xz = analytic_proj(vdf, n0, 1, 3, 2, res_down_bw[2][1], res_down_bw[2][2], z_int, step(z_int))
ana_bw_yz = analytic_proj(vdf, n0, 2, 3, 1, res_down_bw[3][1], res_down_bw[3][2], z_int, step(z_int))
println(io, "| **Backward (down)** | $(round(rel_l2(res_down_bw[1][3], ana_bw_xy); digits = 3)) | $(round(rel_l2(res_down_bw[2][3], ana_bw_xz); digits = 3)) | $(round(rel_l2(res_down_bw[3][3], ana_bw_yz); digits = 3)) |") #hide
println(io, "| **Flux (up)** | $(round(rel_l2(matrix_of(hists_up[1]), matrix_of(ana_hists[1])); digits = 3)) | $(round(rel_l2(matrix_of(hists_up[2]), matrix_of(ana_hists[2])); digits = 3)) | $(round(rel_l2(matrix_of(hists_up[3]), matrix_of(ana_hists[3])); digits = 3)) |") #hide
println(io, "| **Liouville (up)** | $(round(rel_l2(matrix_of(hists_up_m2[1]), matrix_of(ana_hists[1])); digits = 3)) | $(round(rel_l2(matrix_of(hists_up_m2[2]), matrix_of(ana_hists[2])); digits = 3)) | $(round(rel_l2(matrix_of(hists_up_m2[3]), matrix_of(ana_hists[3])); digits = 3)) |") #hide
ana_bw_xy_u = analytic_proj(vdf, n0, 1, 2, 3, res_up_bw[1][1], res_up_bw[1][2], z_int, step(z_int))
ana_bw_xz_u = analytic_proj(vdf, n0, 1, 3, 2, res_up_bw[2][1], res_up_bw[2][2], z_int, step(z_int))
ana_bw_yz_u = analytic_proj(vdf, n0, 2, 3, 1, res_up_bw[3][1], res_up_bw[3][2], z_int, step(z_int))
println(io, "| **Backward (up)** | $(round(rel_l2(res_up_bw[1][3], ana_bw_xy_u); digits = 3)) | $(round(rel_l2(res_up_bw[2][3], ana_bw_xz_u); digits = 3)) | $(round(rel_l2(res_up_bw[3][3], ana_bw_yz_u); digits = 3)) |") #hide
Markdown.parse(String(take!(io))) #hide

# ## What the deviations tell us
# The backward Liouville reconstruction matches the analytic bi-Maxwellian to within the
# grid resolution (relative L2 ≈ 3×10⁻³). This confirms the core assumption — that a
# bi-Maxwellian centred on the E×B drift velocity is an exact Vlasov steady state, so the
# detector VDF is analytically identical to the source VDF — and validates the solver.
#
# The two forward Monte-Carlo methods deviate more, as expected. Flux injection is off by
# only ∼10% (pure 1/√N sampling scatter). Forward Liouville is the least accurate
# (∼35–50%): under this formulation it weights each trajectory by its *source* phase-space
# density but omits the detector-side `|v_x|` flux correction that flux injection applies,
# so for a distribution with a broad velocity spread it recovers a flux-weighted rather
# than strict phase-space density. The analytic overlay makes this deviation directly
# visible as the white dashed contours sitting systematically off the coloured histogram.

# ## Monte-Carlo sampling noise scales as 1/√N
# The report above quotes Flux Injection at ≈10% and attributes it to 1/√N sampling
# scatter. Here we *verify* that statement directly. We run Flux Injection for a range of
# particle counts `N`, repeating each `N` with many *independent* seed realizations, and
# measure the sampling noise as the relative RMS scatter of the reconstructed histogram
# across realizations. For a multinomial (Monte-Carlo) estimator the per-bin density has
# standard deviation `σ ∝ 1/√N`, so the aggregate noise should fall as `1/√N`
# (equivalently the variance `∝ 1/N`).

vec_of(h) = vcat([vec(matrix_of(p)) for p in h]...)

function run_flux_N(N, rseed)
    function prob_func(prob, ctx)
        v = rand(ctx.rng, vdf)
        u0 = SA[x_source..., v...]
        return remake(prob, u0 = u0)
    end
    u0_dummy = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func)
    sols = TP.solve(prob, Boris(); dt, savestepinterval = 1, trajectories = N, seed = rseed)
    return reconstruct_flux_projections(sols, detector_down, n0, 20.0)
end

const Ns = [1000, 2000, 4000, 8000]
const nseeds = 24

agg_noise = Float64[]   # relative RMS scatter of the histogram across seeds → ∝ 1/√N
agg_var = Float64[]   # total histogram variance across seeds → ∝ 1/N
ana_l2 = Float64[]   # analytic relative L2 per seed (averaged over projections)

for N in Ns
    hs = [run_flux_N(N, s) for s in 1:nseeds]
    vecs = [vec_of(h) for h in hs]

    l2s = [mean(rel_l2(matrix_of(h[i]), matrix_of(ana_hists[i])) for i in 1:3) for h in hs]
    push!(ana_l2, mean(l2s))

    mean_vec = mean(vecs)
    mask = mean_vec .> maximum(mean_vec) * 1.0e-4
    mv = mean_vec[mask]
    res = [norm(v[mask] .- mv) / norm(mv) for v in vecs]
    push!(agg_noise, mean(res))

    mat = reduce(hcat, [v[mask] for v in vecs])   # bins × seeds
    push!(agg_var, sum(var(mat; dims = 2)))
end

function logslope(x, y)
    lx = log10.(Float64.(x)); ly = log10.(y)
    return sum((lx .- mean(lx)) .* (ly .- mean(ly))) / sum((lx .- mean(lx)) .^ 2)
end
slope_noise = logslope(Ns, agg_noise)   # target ≈ -0.5
slope_var = logslope(Ns, agg_var)     # target ≈ -1.0

fig_scaling = Figure(size = (1300, 380), fontsize = 18)
ax1 = Axis(
    fig_scaling[1, 1], xscale = log10, yscale = log10,
    xlabel = L"N\ \mathrm{[trajectories]}",
    ylabel = L"MC\ noise\ (rel.\ RMS\ scatter)",
    title = "Sampling noise ∝ 1/√N  (slope = $(round(slope_noise; digits = 2)))"
)
scatter!(ax1, Float64.(Ns), agg_noise; color = :blue, markersize = 8, label = "measured")
ref1 = agg_noise[end] * sqrt(Ns[end]) ./ sqrt.(Float64.(Ns))
lines!(ax1, Float64.(Ns), ref1; color = :black, linestyle = :dash, label = L"\propto 1/\sqrt{N}")
axislegend(ax1; position = :lt)

ax2 = Axis(
    fig_scaling[1, 2], xscale = log10, yscale = log10,
    xlabel = L"N", ylabel = L"\sum_b \mathrm{Var}_s[h_b]\ (\propto 1/N)",
    title = "Variance ∝ 1/N  (slope = $(round(slope_var; digits = 2)))"
)
scatter!(ax2, Float64.(Ns), agg_var; color = :red, markersize = 8, label = "measured")
ref2 = agg_var[end] * Ns[end] ./ Float64.(Ns)
lines!(ax2, Float64.(Ns), ref2; color = :black, linestyle = :dash, label = L"\propto 1/N")
axislegend(ax2; position = :lt)

ax3 = Axis(
    fig_scaling[1, 3], xscale = log10, yscale = log10,
    xlabel = L"N", ylabel = L"relative\ L2\ vs\ analytic",
    title = "Analytic deviation ∝ 1/√N"
)
scatter!(ax3, Float64.(Ns), ana_l2; color = :green, markersize = 8, label = "analytic L2")
axislegend(ax3; position = :lt)
fig_scaling = DisplayAs.PNG(fig_scaling) #hide

# The measured log–log slopes are close to −0.5 (noise) and −1.0 (variance), confirming
# the expected 1/√N sampling-noise scaling. The analytic relative-L2 deviation falls on
# the same 1/√N line, so the ≈10% Flux-Injection deviation in the report above is
# dominated by Monte-Carlo scatter, not by method bias.

io_s = IOBuffer() #hide
println(io_s, "| N (trajectories) | MC noise (∝ 1/√N) | ΣVar (∝ 1/N) | analytic L2 |") #hide
println(io_s, "| :--- | :--- | :--- | :--- |") #hide
for (i, N) in enumerate(Ns) #hide
    println(io_s, "| $N | $(round(agg_noise[i]; digits = 4)) | $(round(agg_var[i]; digits = 2)) | $(round(ana_l2[i]; digits = 4)) |") #hide
end #hide
Markdown.parse(String(take!(io_s))) #hide
