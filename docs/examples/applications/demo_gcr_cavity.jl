# # [Galactic Cosmic Ray Cavity in Earth–Moon Space](@id GCR-Cavity)
#
# This example reproduces the test-particle simulation of the galactic cosmic ray (GCR)
# cavity in Earth–Moon space reported by
# [Shang et al. (2026), *A galactic cosmic ray cavity in Earth-Moon space*](https://doi.org/10.1126/sciadv.adv1908).
#
# Analysis of Chang'E-4 Lunar Lander Neutron and Dosimetry (LND) proton measurements under
# quiet Parker-spiral interplanetary magnetic field (IMF) conditions ($B_x B_y < 0$,
# $B_z \approx 0$ in GSE coordinates) revealed a $\sim 20\%$ reduction in lower-energy GCR
# proton flux ($9$–$34\,\mathrm{MeV}$) in the **prenoon** sector of the lunar orbit
# (lunar phase $\mathrm{LP} = 8\text{–}10\,\mathrm{hM}$), and a shallower reduction at
# higher energies ($42$–$139\,\mathrm{MeV}$). Because a $20\,\mathrm{MeV}$ proton in a
# $5\,\mathrm{nT}$ IMF has a gyroradius $r_g \approx 20\,R_E$ comparable to Earth's
# magnetosphere ($\sim 10\,R_E$), protons spiraling along Parker-spiral IMF lines that
# connect Earth to the prenoon lunar orbit are strongly deflected or magnetically mirrored
# by Earth's magnetic field, carving out a two-sided cavity along the IMF direction.
#
# Rather than forward-tracking millions of particles per $R_E^3$ from upstream and
# downstream reservoirs, we exploit the time-reversibility of the Lorentz force
# ($\mathbf{B} \to -\mathbf{B}$, $\mathbf{v} \to -\mathbf{v}$, $\mathbf{E} \to +\mathbf{E}$)
# with `TestParticle.jl`'s out-of-place `trace` ODE and `EnsembleThreads()` to evaluate
# steady-state connectivity maps and lunar-orbit profiles in seconds. We conclude with a
# direct Liouville-theorem check comparing open-boundary deflection losses against the
# unbounded isotropic limit.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
using TestParticle: c, Rₑ, eV, Proton
import Magnetostatics as MS
using LinearAlgebra: norm, dot, normalize, ×
using Random, Statistics
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## 1. Composite Earth–Moon Electromagnetic Field
#
# We work in Geocentric Solar Ecliptic (**GSE**) coordinates, where $+X$ points toward the
# Sun, $+Y$ points toward dusk, and $+Z$ points toward ecliptic north. The lunar orbit lies
# at $R_{\mathrm{Moon}} = 60\,R_E$ in the $X\text{–}Y$ plane, parameterized by lunar phase
# $\mathrm{LP} \in [0, 24)\,\mathrm{hM}$ starting at full moon (midnight, $-X$):
#
# ```math
# \mathbf{r}_{\mathrm{Moon}}(\mathrm{LP}) = R_{\mathrm{Moon}}
# \left(-\cos\frac{2\pi\,\mathrm{LP}}{24},\; -\sin\frac{2\pi\,\mathrm{LP}}{24},\; 0\right).
# ```
#
# Thus $\mathrm{LP} = 6\,\mathrm{hM}$ is dawn ($-Y$), $\mathrm{LP} = 9\,\mathrm{hM}$ is
# prenoon ($+X, -Y$), $\mathrm{LP} = 12\,\mathrm{hM}$ is noon ($+X$), and
# $\mathrm{LP} = 15\text{–}16\,\mathrm{hM}$ is afternoon ($+X, +Y$).
#
# Following Shang et al. (2026), we superimpose a $5\,\mathrm{nT}$ Parker-spiral IMF
# oriented at $45^\circ$ ($B_x > 0, B_y < 0, B_z = 0$) with Earth's intrinsic dipole field
# (`Magnetostatics.Dipole(TestParticle.BMoment_Earth)`) and the solar wind convection
# electric field $\mathbf{E}_{\mathrm{SW}} = -\mathbf{V}_{\mathrm{SW}} \times \mathbf{B}_{\mathrm{IMF}}$
# for $V_{\mathrm{SW}} = 400\,\mathrm{km\,s^{-1}}$.

const R_MOON = 60.0 * Rₑ
const B_IMF_MAG = 5.0e-9
const B_IMF = SA[B_IMF_MAG / √2, -B_IMF_MAG / √2, 0.0]
const B_HAT = normalize(B_IMF)
const N_PERP = SA[B_HAT[2], -B_HAT[1], 0.0]

const EARTH_DIPOLE = MS.Dipole(TestParticle.BMoment_Earth)
const V_SW = SA[-400.0e3, 0.0, 0.0]
const E_SW = -(V_SW × B_IMF)

@inline B_total(r) = EARTH_DIPOLE(r) + B_IMF
@inline B_rev_earth(r) = -B_total(r)
@inline B_rev_imf(r) = -B_IMF
@inline E_sw_field(r) = E_SW

@inline function lunar_position(lp::Real)
    θ = 2π * Float64(lp) / 24.0
    return SA[-R_MOON * cos(θ), -R_MOON * sin(θ), 0.0]
end

@inline function proton_speed(E_MeV::Real)
    γ = 1.0 + Float64(E_MeV) * 1.0e6 * eV / (Proton.m * c^2)
    return c * sqrt(1.0 - inv(γ^2))
end

@inline function proton_gyroradius_RE(E_MeV::Real)
    γ = 1.0 + Float64(E_MeV) * 1.0e6 * eV / (Proton.m * c^2)
    v = c * sqrt(1.0 - inv(γ^2))
    return γ * Proton.m * v / (Proton.q * B_IMF_MAG * Rₑ)
end

println(
    "20 MeV proton:  v = ", round(proton_speed(20.0) / c; digits = 3),
    " c,  r_g ≈ ", round(proton_gyroradius_RE(20.0); digits = 1), " R_E"
)
println(
    "100 MeV proton: v = ", round(proton_speed(100.0) / c; digits = 3),
    " c,  r_g ≈ ", round(proton_gyroradius_RE(100.0); digits = 1), " R_E"
)

# ### Sample Forward Trajectories: 20 MeV vs. 100 MeV Protons
#
# Before mapping the full cavity, we trace a handful of $20\,\mathrm{MeV}$ and
# $100\,\mathrm{MeV}$ protons launched from the premidnight upstream side along Parker-spiral
# field lines at different impact parameters $d_\perp$ relative to Earth. Protons on field
# lines far from Earth ($d_\perp = 55\,R_E$, passing through the afternoon sector) spiral
# smoothly without perturbation, whereas $20\,\mathrm{MeV}$ protons approaching Earth along
# the prenoon–premidnight flux tube ($d_\perp \lesssim 15\,R_E$) are mirrored or strongly
# scattered by Earth's magnetic field.

param_fwd = prepare(; E = E_sw_field, B = B_total, species = Proton)

function trace_sample_orbits(E_MeV::Float64; d_perps = (-12.0, 0.0, 12.0, 55.0))
    v_mag = proton_speed(E_MeV)
    α = deg2rad(35.0)
    v_par = v_mag * cos(α)
    v_perp = v_mag * sin(α)
    v0 = v_par * B_HAT + v_perp * SA[0.0, 0.0, 1.0]

    cb = TerminateOutside(
        (u, p, t) -> begin
            r2 = u[1]^2 + u[2]^2 + u[3]^2
            r2 < Rₑ^2 || r2 > (150.0 * Rₑ)^2
        end
    )

    sols = Vector{ODESolution}(undef, length(d_perps))
    for (i, dp) in enumerate(d_perps)
        r0 = -105.0 * Rₑ * B_HAT + dp * Rₑ * N_PERP
        prob = ODEProblem(trace, vcat(r0, v0), (0.0, 45.0), param_fwd)
        sols[i] = solve(prob, Tsit5(); callback = cb, reltol = 1.0e-6, abstol = 1.0e-3)
    end
    return sols
end

d_perps_demo = (-12.0, 0.0, 12.0, 55.0)
sols_20 = trace_sample_orbits(20.0; d_perps = d_perps_demo)
sols_100 = trace_sample_orbits(100.0; d_perps = d_perps_demo)

f1 = Figure(size = (1100, 500), fontsize = 14)
θ_circ = range(0, 2π; length = 200)
orbit_x = 60.0 .* cos.(θ_circ)
orbit_y = 60.0 .* sin.(θ_circ)
palette = Makie.wong_colors()

for (col, (E_label, sols)) in enumerate(
        [
            (L"20\ \mathrm{MeV\ Protons}\ (r_g \approx 20\,R_E)", sols_20),
            (L"100\ \mathrm{MeV\ Protons}\ (r_g \approx 47\,R_E)", sols_100),
        ]
    )
    ax = Axis(
        f1[1, col];
        title = E_label,
        xlabel = L"X_{\mathrm{GSE}}\ [R_E]",
        ylabel = col == 1 ? L"Y_{\mathrm{GSE}}\ [R_E]" : "",
        aspect = DataAspect(),
        limits = ((-100.0, 100.0), (-100.0, 100.0)),
    )
    ## Background Parker-spiral IMF lines
    for offset in -120.0:30.0:120.0
        s_line = range(-150.0, 150.0; length = 2)
        lx = [s * B_HAT[1] + offset * N_PERP[1] for s in s_line]
        ly = [s * B_HAT[2] + offset * N_PERP[2] for s in s_line]
        lines!(ax, lx, ly; color = (:gray70, 0.45), linewidth = 1.0)
    end
    lines!(ax, orbit_x, orbit_y; color = :gray30, linestyle = :dash, linewidth = 1.4)
    for (lp_lbl, txt) in [(6, "6 hM"), (9, "9 hM\n(prenoon)"), (12, "12 hM"), (15, "15 hM")]
        p_lp = lunar_position(lp_lbl) ./ Rₑ
        text!(
            ax, p_lp[1] * 1.18, p_lp[2] * 1.18;
            text = txt, align = (:center, :center), fontsize = 11, color = :gray25,
        )
    end
    for (i, sol) in enumerate(sols)
        xs = [u[1] / Rₑ for u in sol.u]
        ys = [u[2] / Rₑ for u in sol.u]
        dp_i = Int(d_perps_demo[i])
        lbl = L"d_\perp = %$dp_i\,R_E"
        lines!(ax, xs, ys; color = palette[i], linewidth = 2.0, label = lbl)
    end
    scatter!(ax, [0.0], [0.0]; color = :dodgerblue, strokecolor = :black, markersize = 10)
    if col == 2
        axislegend(
            ax; position = :lb, framevisible = true, backgroundcolor = (:white, 0.85)
        )
    end
end

f1 = DisplayAs.PNG(f1) #hide

# ## 2. Backward Liouville Tracing of the GCR Cavity
#
# In Shang et al. (2026, Fig. 6), particles are replenished in two source reservoirs
# upstream and downstream along the Parker-spiral IMF and removed upon exiting the
# simulation box or entering Earth's interior ($r < R_E$). By time-reversibility of the
# Lorentz force, a proton is present at $(\mathbf{r}_0, \mathbf{v}_0)$ in steady state if
# and only if its backward trajectory (integrated forward in time with
# $\mathbf{B} \to -\mathbf{B}$ and $\mathbf{E} \to +\mathbf{E}$) reaches one of the two
# source reservoirs ($|s_\parallel| \ge s_{\mathrm{src}}$) within $T_{\mathrm{steady}}$
# without hitting Earth or leaving the transverse box boundaries ($|s_\perp| > L_\perp$ or
# $|z| > L_z$).
#
# Setting the vertical box half-height to $L_z \approx 1.5\,r_g$ ($30\,R_E$ at
# $20\,\mathrm{MeV}$ and $60\,R_E$ at $100\,\mathrm{MeV}$) keeps the control (IMF-only)
# distribution uniform across Earth–Moon space while capturing the transverse deflection
# losses caused by Earth's magnetic obstacle.

const PARAM_REV_EARTH = prepare(; E = E_sw_field, B = B_rev_earth, species = Proton)
const PARAM_REV_IMF = prepare(; E = E_sw_field, B = B_rev_imf, species = Proton)

function sample_connectivity(
        pts::Vector{SVector{3, Float64}}, param, E_MeV::Float64;
        n_dir::Int = 400,
        s_src::Float64 = 160.0 * Rₑ,
        L_perp::Float64 = 240.0 * Rₑ,
        L_z::Float64 = 30.0 * Rₑ,
        T_steady::Float64 = 60.0,
        with_earth::Bool = true,
        seed::Int = 2026,
    )
    v_mag = proton_speed(E_MeV)
    n_pts = length(pts)

    function is_done(u, p, t)
        r = SA[u[1], u[2], u[3]]
        with_earth && dot(r, r) < Rₑ^2 && return true
        abs(dot(r, B_HAT)) >= s_src && return true
        abs(dot(r, N_PERP)) > L_perp && return true
        abs(r[3]) > L_z && return true
        return false
    end

    cb = TerminateOutside(is_done)
    base_prob = ODEProblem(
        trace, SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], (0.0, T_steady), param
    )

    function prob_func(prob, ctx)
        i = ctx.sim_id
        pt_idx = div(i - 1, n_dir) + 1
        r0 = pts[pt_idx]
        rng = Xoshiro(seed + 7919 * i)
        μ = 2.0 * rand(rng) - 1.0
        ϕ = 2π * rand(rng)
        s = sqrt(1.0 - μ^2)
        v0 = SA[v_mag * s * cos(ϕ), v_mag * s * sin(ϕ), v_mag * μ]
        return remake(prob; u0 = vcat(r0, v0))
    end

    function output_func(sol, i)
        u = sol.u[end]
        r = SA[u[1], u[2], u[3]]
        sp = abs(dot(r, B_HAT))
        sn = abs(dot(r, N_PERP))
        az = abs(r[3])
        hit_src = (!with_earth || dot(r, r) >= Rₑ^2) &&
            (sp >= s_src) && (sn <= L_perp) && (az <= L_z)
        return (hit_src, false)
    end

    ens_prob = EnsembleProblem(base_prob; prob_func, output_func, safetycopy = false)
    sols = solve(
        ens_prob, Tsit5(), EnsembleThreads();
        trajectories = n_pts * n_dir,
        callback = cb,
        save_everystep = false,
        save_start = false,
        reltol = 1.0e-5,
        abstol = 1.0e-3,
        maxiters = 20_000,
    )
    hits = reshape(collect(sols.u), n_dir, n_pts)
    return vec(mean(hits; dims = 1))
end

## 2D equatorial GCR maps on X, Y ∈ [-100, 100] R_E
const GRID_XY = range(-100.0, 100.0; length = 29)
grid_pts = vec([SA[x * Rₑ, y * Rₑ, 0.0] for x in GRID_XY, y in GRID_XY])

map_20_raw = reshape(
    sample_connectivity(grid_pts, PARAM_REV_EARTH, 20.0; n_dir = 400, L_z = 30.0 * Rₑ),
    length(GRID_XY), length(GRID_XY),
)
map_100_raw = reshape(
    sample_connectivity(grid_pts, PARAM_REV_EARTH, 100.0; n_dir = 400, L_z = 60.0 * Rₑ),
    length(GRID_XY), length(GRID_XY),
)

## High-statistics lunar-orbit profiles for LP = 8..16 hM
const LPS_SCAN = collect(8.0:1.0:16.0)
lp_pts = [lunar_position(lp) for lp in LPS_SCAN]
const N_DIR_LP = 1600

lp_20_raw = sample_connectivity(
    lp_pts, PARAM_REV_EARTH, 20.0; n_dir = N_DIR_LP, L_z = 30.0 * Rₑ, seed = 3141
)
lp_100_raw = sample_connectivity(
    lp_pts, PARAM_REV_EARTH, 100.0; n_dir = N_DIR_LP, L_z = 60.0 * Rₑ, seed = 2718
)

## Normalize to the unobstructed afternoon value at LP = 16 hM (as in Shang et al. Fig. 5)
norm_20 = lp_20_raw[end]
norm_100 = lp_100_raw[end]
map_20 = map_20_raw ./ norm_20
map_100 = map_100_raw ./ norm_100
lp_20 = lp_20_raw ./ norm_20
lp_100 = lp_100_raw ./ norm_100
err_20 = sqrt.(lp_20_raw .* (1.0 .- lp_20_raw) ./ N_DIR_LP) ./ norm_20
err_100 = sqrt.(lp_100_raw .* (1.0 .- lp_100_raw) ./ N_DIR_LP) ./ norm_100

f2 = Figure(size = (1120, 920), fontsize = 14)

for (row, (E_str, map_val, lp_val, lp_err, col_pt)) in enumerate(
        [
            ("20 MeV Protons", map_20, lp_20, err_20, :royalblue),
            ("100 MeV Protons", map_100, lp_100, err_100, :crimson),
        ]
    )
    ax_map = Axis(
        f2[row, 1];
        title = "$E_str: Equatorial GCR Distribution",
        xlabel = L"X_{\mathrm{GSE}}\ [R_E]",
        ylabel = L"Y_{\mathrm{GSE}}\ [R_E]",
        aspect = DataAspect(),
        limits = ((-100.0, 100.0), (-100.0, 100.0)),
    )
    hm = heatmap!(
        ax_map, GRID_XY, GRID_XY, map_val;
        colormap = :viridis, colorrange = (0.55, 1.08),
    )
    lines!(ax_map, orbit_x, orbit_y; color = :white, linestyle = :dash, linewidth = 1.6)
    ## Highlight the Chang'E-4 LND dayside sector (LP = 7.5 to 16.5 hM)
    lnd_pts = [lunar_position(lp) ./ Rₑ for lp in range(7.5, 16.5; length = 60)]
    lines!(
        ax_map, first.(lnd_pts), [p[2] for p in lnd_pts];
        color = :magenta, linewidth = 2.5,
    )
    scatter!(ax_map, [0.0], [0.0]; color = :black, strokecolor = :white, markersize = 10)
    Colorbar(f2[row, 2], hm; label = "Normalized Count")

    ax_lp = Axis(
        f2[row, 3];
        title = "$E_str: Lunar Orbit Profile",
        xlabel = "Lunar Phase [hM]",
        ylabel = "Normalized Count (rel. to 16 hM)",
        limits = ((7.5, 16.5), (0.65, 1.12)),
        xticks = 8:2:16,
    )
    vspan!(ax_lp, 8.0, 10.0; color = (:gray, 0.18), label = "Prenoon cavity (8–10 hM)")
    hlines!(ax_lp, [1.0]; color = :gray45, linestyle = :dash, linewidth = 1.2)
    errorbars!(ax_lp, LPS_SCAN, lp_val, lp_err; color = col_pt, whiskerwidth = 6)
    scatterlines!(
        ax_lp, LPS_SCAN, lp_val;
        color = col_pt, marker = :utriangle, markersize = 11, linewidth = 2.0,
        label = "TestParticle.jl",
    )
    axislegend(ax_lp; position = :rb, framevisible = false)
end

colsize!(f2.layout, 1, Relative(0.44))
colsize!(f2.layout, 3, Relative(0.46))

f2 = DisplayAs.PNG(f2) #hide

# As in Fig. 5 of Shang et al. (2026), a two-sided GCR cavity extends along the
# Parker-spiral IMF direction into both the prenoon ($\mathrm{LP} = 8\text{–}10\,\mathrm{hM}$)
# and premidnight ($\mathrm{LP} = 20\text{–}22\,\mathrm{hM}$) sectors. At $20\,\mathrm{MeV}$
# the normalized count drops by $\sim 20\%$ near $\mathrm{LP} = 9\,\mathrm{hM}$ and recovers
# to unity in the afternoon ($\mathrm{LP} = 13\text{–}16\,\mathrm{hM}$), while the
# $100\,\mathrm{MeV}$ cavity is noticeably shallower because of the larger proton gyroradius.

# ## 3. IMF-Angle Dependence and Liouville-Theorem Check
#
# Finally, we examine two deeper physical aspects of the cavity:
#
# 1. **Unbounded magnetospheric mirroring vs. IMF angle $\varphi$**: Along the lunar orbit,
#    let $\varphi \in [0^\circ, 180^\circ]$ be the angle between $\mathbf{B}_{\mathrm{IMF}}$
#    and the Earth–Moon vector $\mathbf{r}_{\mathrm{Moon}}$. Tracing backward in an
#    *unbounded* domain ($r_{\max} = 200\,R_E$, no vertical wall), we measure the fraction
#    $1 - P_{\mathrm{mirrored}}$ of trajectories that arrive unobstructed without being
#    reflected ($v_\parallel$ sign reversal inside $r < 25\,R_E$) by Earth's magnetic field.
#    This quantity is independent of box boundaries and reproduces the concave $\varphi$
#    dependence of Shang et al. (2026, Fig. 4): minimum near $\varphi = 0^\circ$ and
#    $180^\circ$ (field-aligned Earth–Moon connection) and maximum near $\varphi = 90^\circ$.
# 2. **Liouville-theorem boundary check**: In a static electromagnetic field with a uniform,
#    isotropic, *unbounded* bath at infinity, Liouville's theorem requires $f(\mathbf{r},
#    \mathbf{v}) = f_0$ on every trajectory that escapes to infinity, even if it is mirrored
#    or deflected by Earth's magnetosphere. Because Earth's $1\,R_E$ solid body absorbs
#    $< 0.2\%$ of $20\,\mathrm{MeV}$ protons at $60\,R_E$, increasing the vertical box
#    half-height from $L_z = 30\,R_E$ to $L_z = 80\,R_E$ (with $T_{\mathrm{steady}} = 180\,\mathrm{s}$)
#    causes the omnidirectional cavity to vanish — demonstrating that in collisionless
#    test-particle models, a steady-state omnidirectional density drop at $60\,R_E$ requires
#    either open transverse boundaries, an absorbing magnetospheric sink, or a directional
#    (pitch-angle / field-of-view) anisotropy.

function sample_unobstructed_vs_phi(E_MeV::Float64; n_dir::Int = 1000, seed::Int = 9001)
    v_mag = proton_speed(E_MeV)
    ## LP = 9..21 hM sweeps monotonically from φ = 0° to φ = 180°
    lps_half = collect(9.0:1.0:21.0)
    pts = [lunar_position(lp) for lp in lps_half]
    phis = [
        rad2deg(acos(clamp(dot(normalize(p), B_HAT), -1.0, 1.0))) for p in pts
    ]

    cb = TerminateOutside(
        (u, p, t) -> begin
            r2 = u[1]^2 + u[2]^2 + u[3]^2
            r2 < Rₑ^2 || r2 > (200.0 * Rₑ)^2
        end
    )
    base_prob = ODEProblem(
        trace, SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], (0.0, 60.0), PARAM_REV_EARTH
    )

    function prob_func(prob, ctx)
        i = ctx.sim_id
        pt_idx = div(i - 1, n_dir) + 1
        r0 = pts[pt_idx]
        rng = Xoshiro(seed + 7919 * i)
        μ = 2.0 * rand(rng) - 1.0
        ϕ = 2π * rand(rng)
        s = sqrt(1.0 - μ^2)
        v0 = SA[v_mag * s * cos(ϕ), v_mag * s * sin(ϕ), v_mag * μ]
        return remake(prob; u0 = vcat(r0, v0))
    end

    function output_func(sol, i)
        v0 = SA[sol.u[1][4], sol.u[1][5], sol.u[1][6]]
        vend = SA[sol.u[end][4], sol.u[end][5], sol.u[end][6]]
        rmin2 = minimum(u[1]^2 + u[2]^2 + u[3]^2 for u in sol.u)
        mirrored = (dot(v0, B_HAT) * dot(vend, B_HAT) < 0.0) && (rmin2 < (25.0 * Rₑ)^2)
        hit_earth = rmin2 < Rₑ^2
        return (!(mirrored || hit_earth), false)
    end

    ens_prob = EnsembleProblem(base_prob; prob_func, output_func, safetycopy = false)
    sols = solve(
        ens_prob, Tsit5(), EnsembleThreads();
        trajectories = length(pts) * n_dir,
        callback = cb,
        save_everystep = true,
        reltol = 1.0e-5,
        abstol = 1.0e-3,
        maxiters = 20_000,
    )
    unobstructed = vec(mean(reshape(collect(sols.u), n_dir, length(pts)); dims = 1))
    ord = sortperm(phis)
    return phis[ord], unobstructed[ord]
end

phi_20, trans_20 = sample_unobstructed_vs_phi(20.0; n_dir = 800, seed = 9001)
phi_100, trans_100 = sample_unobstructed_vs_phi(100.0; n_dir = 800, seed = 9002)

lp_20_unbounded = sample_connectivity(
    lp_pts, PARAM_REV_EARTH, 20.0;
    n_dir = N_DIR_LP, L_z = 80.0 * Rₑ, T_steady = 180.0, seed = 4242,
)
lp_20_unbounded_norm = lp_20_unbounded ./ lp_20_unbounded[end]

f3 = Figure(size = (1120, 440), fontsize = 14)

ax_phi = Axis(
    f3[1, 1];
    title = "Unobstructed Fraction vs. IMF Angle φ (Unbounded)",
    xlabel = L"\varphi = \angle(\mathbf{B}_{\mathrm{IMF}}, \mathbf{r}_{\mathrm{Moon}})" *
        L"\ [^\circ]",
    ylabel = L"1 - P_{\mathrm{mirrored}}",
    limits = ((-5.0, 185.0), (0.75, 1.03)),
    xticks = 0:30:180,
)
scatterlines!(
    ax_phi, phi_20, trans_20;
    color = :royalblue, linewidth = 2.0, markersize = 8, label = "20 MeV",
)
scatterlines!(
    ax_phi, phi_100, trans_100;
    color = :crimson, linewidth = 2.0, markersize = 8, label = "100 MeV",
)
axislegend(ax_phi; position = :cb, framevisible = false)

ax_liouv = Axis(
    f3[1, 2];
    title = "20 MeV Liouville Check: Open Box vs. Unbounded Limit",
    xlabel = "Lunar Phase [hM]",
    ylabel = "Normalized Count (rel. to 16 hM)",
    limits = ((7.5, 16.5), (0.70, 1.12)),
    xticks = 8:2:16,
)
vspan!(ax_liouv, 8.0, 10.0; color = (:gray, 0.18))
hlines!(ax_liouv, [1.0]; color = :gray45, linestyle = :dash, linewidth = 1.2)
scatterlines!(
    ax_liouv, LPS_SCAN, lp_20;
    color = :royalblue, marker = :utriangle, markersize = 10, linewidth = 2.0,
    label = L"Open box ($L_z = 30\,R_E$, $T = 60\,\mathrm{s}$)",
)
scatterlines!(
    ax_liouv, LPS_SCAN, lp_20_unbounded_norm;
    color = :seagreen, marker = :circle, markersize = 9, linewidth = 2.0,
    label = L"Unbounded limit ($L_z = 80\,R_E$, $T = 180\,\mathrm{s}$)",
)
axislegend(ax_liouv; position = :rb, framevisible = false)

f3 = DisplayAs.PNG(f3) #hide
