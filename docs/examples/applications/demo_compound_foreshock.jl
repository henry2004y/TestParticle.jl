# # Compound Electron Acceleration at Planetary Foreshocks
#
# This example demonstrates **compound electron acceleration** at a quasi-perpendicular
# planetary bow shock, inspired by
# [Shi et al. (2025), *Compound electron acceleration at planetary foreshocks*](https://www.nature.com/articles/s41467-024-55464-8).
# Solar wind electrons have typical thermal energies of $10$–$100\,\mathrm{eV}$, yet
# spacecraft observations at Earth's bow shock frequently reveal power-law electron tails
# extending above $100\,\mathrm{keV}$.
#
# No single mechanism produces this multi-order-of-magnitude energization on its own:
#
# 1. **Shock Drift Acceleration (SDA)** at the quasi-perpendicular bow shock
#    ($\theta_{Bn} \approx 80^\circ$) energizes reflecting electrons, but without
#    scattering they quickly enter the magnetic loss cone and escape after a few bounces.
# 2. **Multi-mode resonant pitch-angle scattering** localized around the compressional
#    boundaries—electrostatic turbulence (EST, $\lesssim 1\,\mathrm{keV}$),
#    high-frequency whistler-mode waves (WW, $\sim 0.2$–$10\,\mathrm{keV}$), and
#    low-frequency magnetosonic waves (MSW, $\gtrsim 2\,\mathrm{keV}$)—scatters electrons
#    out of the loss cone so they remain trapped and undergo repeated shock reflections.
# 3. **Fermi and betatron acceleration** from a train of approaching foreshock transient
#    compressional boundaries ($n_t = 3$) and the resulting core magnetic compression
#    further boost trapped electrons past $100\,\mathrm{keV}$ ($\gamma > 1.2$).
#
# Whereas Shi et al. (2025) modeled these effects with a 0D bounce-averaged Monte-Carlo
# map, here we build a self-consistent **relativistic guiding-center** model in
# `TestParticle.jl` (`RelativisticGuidingCenter`) using Maxwell-consistent analytical
# $\mathbf{E}(\mathbf{x}, t)$ and $\mathbf{B}(\mathbf{x}, t)$ fields coupled with an
# `OrdinaryDiffEq` `DiscreteCallback` for boundary-localized pitch-angle diffusion.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
import TestParticle as TP
using TestParticle: Rₑ, eV, c
using PlasmaPopulations: BiMaxwellian
using LinearAlgebra: norm, tr
using Random, Statistics, Printf, Markdown
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Maxwell-Consistent Multi-Transient Foreshock Fields
#
# We model the foreshock region along the shock-normal direction $x$ with a stationary
# quasi-perpendicular bow shock at $x = 0$ ($\theta_{Bn} = 80^\circ$, $R_{\mathrm{bs}} = 4$)
# and a train of $n_t$ approaching foreshock transient compressional boundaries at
# $x_k(t) = k\,\Delta L - v_f t$ ($k = 1, \dots, n_t$, $R_{\mathrm{ft}} = 4$). Keeping
# $B_x = B_0 \cos\theta_{Bn}$ constant guarantees $\nabla\cdot\mathbf{B} = 0$ identically:
#
# ```math
# \mathbf{B}(x, t) = \begin{bmatrix} B_{x0} \\ 0 \\ B_z(x, t) \end{bmatrix},
# \qquad
# B_z(x, t) = B_{z0}\,(1 + \beta t)\left[
#   1 + (R_{\mathrm{bs}} - 1)\,\sigma\!\left(\frac{-x}{w_{\mathrm{bs}}}\right)
#     + (R_{\mathrm{ft}} - 1)\sum_{k=1}^{n_t}
#       \operatorname{sech}^2\!\left(\frac{x - x_k(t)}{w_{\mathrm{ft}}}\right)
# \right],
# ```
#
# where $\sigma(s) = \frac{1}{2}(1 + \tanh s)$ and $1 + \beta t$ represents the core
# magnetic compression (reaching $3\times$ over $t_{\mathrm{end}}$). Faraday's law
# $\nabla\times\mathbf{E} = -\partial\mathbf{B}/\partial t$ reduces to
# $\partial E_y / \partial x = -\partial B_z / \partial t$, which integrates analytically in
# closed form on top of the steady solar-wind convection electric field
# $E_{y0} = -v_{\mathrm{sw}} B_{z0}$.

const B₀ = 5.0e-9              # [T] ambient upstream magnetic field
const θ_Bn = deg2rad(80.0)     # [rad] MMS quasi-perpendicular bow shock angle
const Bx₀ = B₀ * cos(θ_Bn)     # [T] constant normal component
const Bz₀ = B₀ * sin(θ_Bn)     # [T] ambient tangential component
const v_sw = 650.0e3           # [m/s] MMS shock-normal flow speed
const Ey₀ = -v_sw * Bz₀        # [V/m] steady convection electric field

const R_bs = 4.0               # bow shock compression ratio
const R_ft = 4.0               # foreshock transient boundary compression ratio
const ΔL = 1.8 * Rₑ            # [m] spacing between consecutive transients
const v_f = 280.0e3            # [m/s] foreshock transient approach speed
const w_bs = 0.05 * Rₑ         # [m] bow shock ramp width
const w_ft = 0.08 * Rₑ         # [m] foreshock transient boundary half-width
const w_wave = 0.25 * Rₑ       # [m] boundary wave-scattering layer half-width
const t_end = 115.0            # [s] simulation duration (~3 transient collapses)
const β_comp = 2.0 / t_end     # [1/s] core compression rate (3× B over t_end)

@inline σ(s) = 0.5 * (1.0 + tanh(s))
@inline log2cosh(s) = abs(s) + log1p(exp(-2.0 * abs(s)))

function make_foreshock_fields(; n_t::Int = 3, β = 0.0, vf = 0.0)
    L_init = ntuple(k -> k * ΔL, n_t)
    function B(x, t)
        s_bs = σ(-x[1] / w_bs)
        bumps = 0.0
        for k in 1:n_t
            xk = L_init[k] - vf * t
            bumps += sech((x[1] - xk) / w_ft)^2
        end
        Bz = Bz₀ * (1.0 + β * t) * (1.0 + (R_bs - 1.0) * s_bs + (R_ft - 1.0) * bumps)
        return SA[Bx₀, 0.0, Bz]
    end
    function E(x, t)
        I_bs = 0.5 * (x[1] - w_bs * (log2cosh(-x[1] / w_bs) - log2cosh(0.0)))
        I_bumps = 0.0
        d_bumps = 0.0
        for k in 1:n_t
            xk = L_init[k] - vf * t
            I_bumps += w_ft * (tanh((x[1] - xk) / w_ft) - tanh(-xk / w_ft))
            d_bumps += sech((x[1] - xk) / w_ft)^2 - sech(-xk / w_ft)^2
        end
        Ey_beta = -Bz₀ * β * (x[1] + (R_bs - 1.0) * I_bs + (R_ft - 1.0) * I_bumps)
        Ey_ft = -Bz₀ * (1.0 + β * t) * (R_ft - 1.0) * vf * d_bumps
        return SA[0.0, Ey₀ + Ey_beta + Ey_ft, 0.0]
    end
    return B, E, L_init
end

# Verify that the multi-transient field pair satisfies Maxwell's equations everywhere:

let
    Bfunc, Efunc, _ = make_foreshock_fields(; n_t = 3, β = β_comp, vf = v_f)
    Bf, Ef = TP.Field(Bfunc), TP.Field(Efunc)
    t = 30.0
    for x_re in (0.0, 1.0, 2.5, 4.0)
        x = SA[x_re * Rₑ, 0.0, 0.0]
        JB = TP.jacobian(Bf, x, t)
        JE = TP.jacobian(Ef, x, t)
        curlE = SA[JE[3, 2] - JE[2, 3], JE[1, 3] - JE[3, 1], JE[2, 1] - JE[1, 2]]
        dBdt = TP.derivative_t(Bf, x, t)
        @printf(
            "x = %.1f Rₑ: |∇⋅B| = %.2e, ‖∇×E + ∂B/∂t‖ / ‖∂B/∂t‖ = %.2e\n",
            x_re, abs(tr(JB)), norm(curlE + dBdt) / norm(dBdt)
        )
    end
end

# ## Boundary-Localized Multi-Mode Pitch-Angle Diffusion
#
# Following Fig. 3 of Shi et al. (2025), wave activity is concentrated near the bow shock
# and foreshock transient compressional boundaries, and spans three overlapping energy
# bands:
#
# - **EST** (electrostatic turbulence): active below $\sim 1\,\mathrm{keV}$,
# - **WW** (high-frequency whistler-mode waves): bridges $\sim 0.2$–$10\,\mathrm{keV}$,
# - **MSW** (low-frequency magnetosonic waves): scatters $\gtrsim 2\,\mathrm{keV}$ up to
#   relativistic energies.

@inline function wave_localization(x1, t, L_init, vf)
    w = sech(x1 / w_wave)^2
    for xk0 in L_init
        xk = xk0 - vf * t
        w += sech((x1 - xk) / w_wave)^2
    end
    return min(w, 1.0)
end

@inline function D_αα(E_eV; est = false, ww = false, msw = false)
    logE = log10(max(E_eV, 1.0))
    D = 0.0
    if est
        D += 2.0 * exp(-((logE - 2.0) / 0.7)^2)
    end
    if ww
        D += 1.6 * exp(-((logE - 3.3) / 0.65)^2)
    end
    if msw
        D += 1.4 * exp(-((logE - 4.7) / 0.9)^2)
    end
    return D
end

# Visualize the compressing multi-transient magnetic profile and the three wave bands:

f1 = Figure(size = (1050, 420), fontsize = 14)

ax_B = Axis(
    f1[1, 1]; title = "Multi-transient foreshock magnetic field (n_t = 3)",
    xlabel = L"x\ [R_E]", ylabel = L"|\mathbf{B}|\ [\mathrm{nT}]",
)
let
    Bfunc, _, _ = make_foreshock_fields(; n_t = 3, β = β_comp, vf = v_f)
    xs = range(-0.4, 6.2; length = 600)
    for (t_val, col, lbl) in (
            (0.0, :steelblue, "t = 0 s"),
            (55.0, :darkorange, "t = 55 s"),
            (t_end, :crimson, "t = $(Int(t_end)) s"),
        )
        Bvals = [norm(Bfunc(SA[x * Rₑ, 0.0, 0.0], t_val)) * 1.0e9 for x in xs]
        lines!(ax_B, xs, Bvals; color = col, linewidth = 2, label = lbl)
    end
    vlines!(ax_B, [0.0]; color = :gray40, linestyle = :dash, label = "Bow shock")
    axislegend(ax_B; position = :rt)
end

ax_D = Axis(
    f1[1, 2]; title = "Boundary wave pitch-angle diffusion bands",
    xlabel = L"E\ [\mathrm{keV}]", ylabel = L"D_{\alpha\alpha}\ [\mathrm{rad}^2/\mathrm{s}]",
    xscale = log10, yscale = log10, limits = ((0.01, 300.0), (1.0e-3, 5.0)),
)
let
    Es_eV = 10.0 .^ range(1.0, 5.5; length = 300)
    Es_keV = Es_eV ./ 1.0e3
    lines!(
        ax_D, Es_keV, [D_αα(E; est = true) for E in Es_eV];
        color = :crimson, linewidth = 2, label = "EST"
    )
    lines!(
        ax_D, Es_keV, [D_αα(E; ww = true) for E in Es_eV];
        color = :forestgreen, linewidth = 2, label = "WW"
    )
    lines!(
        ax_D, Es_keV, [D_αα(E; msw = true) for E in Es_eV];
        color = :royalblue, linewidth = 2, label = "MSW"
    )
    lines!(
        ax_D, Es_keV, [D_αα(E; est = true, ww = true, msw = true) for E in Es_eV];
        color = :black, linestyle = :dash, linewidth = 2, label = "Total"
    )
    axislegend(ax_D; position = :rb)
end

f1 = DisplayAs.PNG(f1) #hide

# ## Relativistic Guiding-Center Tracing with Pitch-Angle Scattering
#
# In `RelativisticGuidingCenter`, the ODE state is
# $\mathbf{y} = (X_x, X_y, X_z, u_\parallel, \mu)$ with parallel proper velocity
# $u_\parallel = \gamma v_\parallel$ and relativistic magnetic moment
# $\mu = p_\perp^2 / (2 m_e B) = m_e u_\perp^2 / (2 B)$, so a scattering kick resets both
# $u_\parallel$ and $\mu$ in the state. Each trajectory carries its own `Xoshiro` RNG in
# its parameters so `EnsembleThreads()` remains thread-safe and reproducible.
#
# Isotropic pitch-angle diffusion on the velocity sphere obeys the singularity-free Itô SDE
# in $\mu_\alpha = \cos\alpha \in [-1, 1]$:
#
# ```math
# d\mu_\alpha = -2 D_{\alpha\alpha}(x, t, E)\,\mu_\alpha\,dt
#   + \sqrt{2 D_{\alpha\alpha}(x, t, E)\,(1 - \mu_\alpha^2)\,dt}\;\mathcal{N}(0, 1),
# ```
#
# applied at the end of each adaptive ODE step via a stateless `DiscreteCallback` with
# reflecting boundaries at $\mu_\alpha = \pm 1$, conserving total kinetic energy
# $K = m_e c^2(\gamma - 1)$ across each elastic scattering kick.

@inline function gc_energy_eV(u_par, μ, Bmag)
    m = Electron.m
    u² = u_par^2 + 2.0 * μ * Bmag / m
    γ = sqrt(1.0 + u² / c^2)
    return m * c^2 * (γ - 1.0) / eV
end

function make_scatter_callback(L_init, vf; est = false, ww = false, msw = false)
    m = Electron.m
    condition(u, t, integrator) = t > integrator.tprev
    function affect!(integrator)
        dt = integrator.t - integrator.tprev
        dt <= 0.0 && return
        u = integrator.u
        t = integrator.t
        wloc = wave_localization(u[1], t, L_init, vf)
        wloc < 1.0e-3 && return
        p = integrator.p
        Bmag = norm(p.B(SA[u[1], u[2], u[3]], t))
        u_par = u[4]
        u_perp² = max(2.0 * u[5] * Bmag / m, 0.0)
        u_tot = sqrt(u_par^2 + u_perp²)
        u_tot == 0.0 && return
        E_eV = m * c^2 * (sqrt(1.0 + u_tot^2 / c^2) - 1.0) / eV
        D = wloc * D_αα(E_eV; est, ww, msw)
        D <= 0.0 && return
        μ_α = clamp(u_par / u_tot, -1.0, 1.0)
        sin²α = max(1.0 - μ_α^2, 0.0)
        μ_new = μ_α - 2.0 * D * μ_α * dt + sqrt(2.0 * D * sin²α * dt) * randn(p.rng)
        if μ_new > 1.0
            μ_new = 2.0 - μ_new
        elseif μ_new < -1.0
            μ_new = -2.0 - μ_new
        end
        μ_new = clamp(μ_new, -1.0, 1.0)
        u_par_new = u_tot * μ_new
        u_perp²_new = u_tot^2 * max(1.0 - μ_new^2, 0.0)
        μ_mag_new = 0.5 * m * u_perp²_new / Bmag
        integrator.u = SVector{5}(u[1], u[2], u[3], u_par_new, μ_mag_new)
        return
    end
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

# ## Single-Electron Compound Trajectory
#
# Before running the full ensemble study, let us inspect a single electron undergoing
# compound acceleration across $n_t = 3$ approaching foreshock transients. We sample a
# small batch of $50\,\mathrm{eV}$ solar wind electrons, pick the one with the largest
# energy gain, and record $(x, K, \alpha, |\mathbf{B}|)$ along its orbit.

function trace_recorded_electron(
        u0_gc, p_init, L_init, vf; est = true, ww = true, msw = true
    )
    m = Electron.m
    ts = Float64[0.0]
    xs = Float64[u0_gc[1] / Rₑ]
    B0 = norm(p_init.B(SA[u0_gc[1], u0_gc[2], u0_gc[3]], 0.0))
    Ks = Float64[gc_energy_eV(u0_gc[4], u0_gc[5], B0) / 1.0e3]
    u_perp0 = sqrt(max(2.0 * u0_gc[5] * B0 / m, 0.0))
    αs = Float64[rad2deg(atan(u_perp0, u0_gc[4]))]
    Bs = Float64[B0 * 1.0e9]

    p = (; p_init..., rng = copy(p_init.rng))
    cb_scat = make_scatter_callback(L_init, vf; est, ww, msw)
    function record_affect!(integrator)
        u = integrator.u
        t = integrator.t
        Bmag = norm(integrator.p.B(SA[u[1], u[2], u[3]], t))
        u_perp = sqrt(max(2.0 * u[5] * Bmag / m, 0.0))
        push!(ts, t)
        push!(xs, u[1] / Rₑ)
        push!(Ks, gc_energy_eV(u[4], u[5], Bmag) / 1.0e3)
        push!(αs, rad2deg(atan(u_perp, u[4])))
        push!(Bs, Bmag * 1.0e9)
        return
    end
    cb_rec = DiscreteCallback(
        (u, t, integrator) -> t > integrator.tprev, record_affect!;
        save_positions = (false, false)
    )
    L_outer = L_init[end]
    isoutside(u, p, t) = u[1] < -4.0 * w_bs || u[1] > (L_outer - vf * t) + 4.0 * w_ft

    prob = TraceProblem{RelativisticGuidingCenter}(u0_gc, (0.0, t_end), p)
    solve(
        prob, Tsit5();
        isoutside, callback = CallbackSet(cb_scat, cb_rec),
        save_everystep = false, dtmax = 0.15, reltol = 1.0e-5, abstol = 1.0e-5,
    )
    return (; ts, xs, Ks, αs, Bs)
end

best_traj, L_init_3 = let
    Bfunc, Efunc, L_init = make_foreshock_fields(; n_t = 3, β = β_comp, vf = v_f)
    p_field = prepare(; E = Efunc, B = Bfunc, species = Electron)
    rng = Xoshiro(2025)
    best_E = 0.0
    best_data = nothing
    for _ in 1:250
        x0 = SA[(0.15 + 0.75 * rand(rng)) * ΔL, 0.0, 0.0]
        v0 = rand(rng, BiMaxwellian(Electron; T = 50.0))
        γ0 = inv(sqrt(1.0 - sum(abs2, v0) / c^2))
        xv0 = SA[x0[1], x0[2], x0[3], γ0 * v0[1], γ0 * v0[2], γ0 * v0[3]]
        sgc = TP.full_to_gc(xv0, p_field; model = RelativisticGuidingCenter())
        p_init = (; p_field..., rng = Xoshiro(rand(rng, UInt64)))
        tr = trace_recorded_electron(sgc, p_init, L_init, v_f)
        if tr.Ks[end] > best_E
            best_E = tr.Ks[end]
            best_data = tr
        end
    end
    @printf(
        "Representative electron: %.1f eV → %.2f keV (×%.0f energy gain)\n",
        best_data.Ks[1] * 1.0e3, best_data.Ks[end], best_data.Ks[end] / best_data.Ks[1]
    )
    best_data, L_init
end

f2 = Figure(size = (950, 740), fontsize = 14)

ax_x = Axis(
    f2[1, 1]; ylabel = L"x\ [R_E]", limits = (nothing, (-0.3, 5.8)),
    title = "Compound acceleration of a single electron (n_t = 3)",
)
lines!(ax_x, best_traj.ts, best_traj.xs; color = :steelblue, linewidth = 1.6, label = "Electron")
for (k, xk0) in enumerate(L_init_3)
    x_tr = [(xk0 - v_f * t) / Rₑ for t in best_traj.ts]
    x_tr = [v >= -0.2 ? v : NaN for v in x_tr]
    lines!(
        ax_x, best_traj.ts, x_tr;
        color = :crimson, linestyle = :dash, linewidth = 1.3,
        label = k == 1 ? "Foreshock transients" : nothing,
    )
end
hlines!(ax_x, [0.0]; color = :gray30, linestyle = :dot, linewidth = 1.5, label = "Bow shock")
axislegend(ax_x; position = :rt, orientation = :horizontal)

ax_K = Axis(f2[2, 1]; ylabel = L"K\ [\mathrm{keV}]", yscale = log10)
lines!(ax_K, best_traj.ts, best_traj.Ks; color = :darkgreen, linewidth = 1.8)

ax_α = Axis(f2[3, 1]; ylabel = L"\alpha\ [\mathrm{deg}]", limits = (nothing, (0, 180)))
α_lc = rad2deg(asin(sqrt(1.0 / R_bs)))
hspan!(ax_α, 0.0, α_lc; color = (:crimson, 0.12))
hspan!(ax_α, 180.0 - α_lc, 180.0; color = (:crimson, 0.12))
lines!(ax_α, best_traj.ts, best_traj.αs; color = :purple, linewidth = 1.0)

ax_Bm = Axis(f2[4, 1]; xlabel = L"t\ [\mathrm{s}]", ylabel = L"|\mathbf{B}|\ [\mathrm{nT}]")
lines!(ax_Bm, best_traj.ts, best_traj.Bs; color = :darkorange, linewidth = 1.2)

linkxaxes!(ax_x, ax_K, ax_α, ax_Bm)
hidexdecorations!.((ax_x, ax_K, ax_α); grid = false)

f2 = DisplayAs.PNG(f2) #hide

# As the first transient collapses into the bow shock, the electron is caught by the
# subsequent approaching transients; boundary-localized wave scattering keeps its pitch
# angle outside the loss cones (red bands, $\alpha < 30^\circ$ and $\alpha > 150^\circ$),
# enabling repeated SDA/Fermi jumps on top of steady betatron heating.

# ## Ensemble Ablation and Multi-Transient Study
#
# Next we reproduce both the progressive mechanism comparison (Fig. 4 of Shi et al., 2025)
# and the sensitivity to the number of foreshock transients $n_t \in \{1, 2, 3\}$ (Fig. 8a)
# by tracing $10{,}000$ electrons per stage from a $50\,\mathrm{eV}$ solar wind Maxwellian.

function run_ensemble(;
        n_t::Int = 3, est = false, ww = false, msw = false, β = 0.0, vf = 0.0,
        trajectories = 10_000, seed = 1234,
    )
    Bfunc, Efunc, L_init = make_foreshock_fields(; n_t, β, vf)
    tspan = (0.0, t_end)

    u0_dummy = SA[0.5 * ΔL, 0.0, 0.0, 1.0e6, 1.0e6, 0.0]
    p_field = prepare(; E = Efunc, B = Bfunc, species = Electron)
    p0 = (; p_field..., rng = Xoshiro(seed))

    L_outer = L_init[end]
    isoutside(u, p, t) = u[1] < -4.0 * w_bs || u[1] > (L_outer - vf * t) + 4.0 * w_ft

    function prob_func(prob, ctx)
        rng = ctx.rng
        x0 = SA[(0.15 + 0.75 * rand(rng)) * ΔL, 0.0, 0.0]
        v0 = rand(rng, BiMaxwellian(Electron; T = 50.0))
        γ0 = inv(sqrt(1.0 - sum(abs2, v0) / c^2))
        xv0 = SA[x0[1], x0[2], x0[3], γ0 * v0[1], γ0 * v0[2], γ0 * v0[3]]
        ## The problem converts the full-orbit state to its guiding center.
        return remake(prob; u0 = xv0, p = (; p_field..., rng = Xoshiro(rand(rng, UInt64))))
    end

    cb_scat = (est || ww || msw) ?
        make_scatter_callback(L_init, vf; est, ww, msw) : nothing
    prob = TraceProblem{RelativisticGuidingCenter}(u0_dummy, tspan, p0; prob_func)
    sols = solve(
        prob, Tsit5(), EnsembleThreads();
        isoutside, callback = cb_scat, trajectories, seed,
        save_everystep = false, dtmax = 0.15, reltol = 1.0e-5, abstol = 1.0e-5,
    )

    E_final = Vector{Float64}(undef, trajectories)
    for (k, sol) in enumerate(sols.u)
        u_end = sol.u[end]
        t_k = sol.t[end]
        B_end = norm(sol.prob.p.B(SA[u_end[1], u_end[2], u_end[3]], t_k))
        E_final[k] = gc_energy_eV(u_end[4], u_end[5], B_end)
    end
    return E_final
end

const N_particles = 10_000

## Initial solar wind energy distribution (50 eV Maxwellian)
E_sw = let
    rng = Xoshiro(1234)
    [
        let v = rand(rng, BiMaxwellian(Electron; T = 50.0))
            0.5 * Electron.m * sum(abs2, v) / eV
        end
            for _ in 1:N_particles
    ]
end

E_full_3 = run_ensemble(;
    n_t = 3, est = true, ww = true, msw = true, β = β_comp, vf = v_f,
    trajectories = N_particles,
)

cases = [
    ("SDA only", run_ensemble(; trajectories = N_particles)),
    ("SDA + EST", run_ensemble(; est = true, trajectories = N_particles)),
    ("SDA + EST + WW", run_ensemble(; est = true, ww = true, trajectories = N_particles)),
    (
        "SDA + EST + WW + MSW",
        run_ensemble(; est = true, ww = true, msw = true, trajectories = N_particles),
    ),
    ("Full (+ Fermi + Betatron)", E_full_3),
]

nt_cases = [
    (
        "n_t = 1",
        run_ensemble(;
            n_t = 1, est = true, ww = true, msw = true, β = β_comp, vf = v_f,
            trajectories = N_particles,
        ),
    ),
    (
        "n_t = 2",
        run_ensemble(;
            n_t = 2, est = true, ww = true, msw = true, β = β_comp, vf = v_f,
            trajectories = N_particles,
        ),
    ),
    ("n_t = 3", E_full_3),
]

io = IOBuffer() #hide
println(io, "| Stage | Median [eV] | 95th pct [keV] | 99th pct [keV] | Max [keV] | Max γ |") #hide
println(io, "| ----- | ----------- | -------------- | -------------- | --------- | ----- |") #hide
for (label, Es) in cases #hide
    emax = maximum(Es) #hide
    γmax = 1.0 + emax * eV / (Electron.m * c^2) #hide
    @printf( #hide
        io, "| %s | %.0f | %.2f | %.2f | %.1f | %.3f |\n", #hide
        label, median(Es), quantile(Es, 0.95) / 1.0e3, #hide
        quantile(Es, 0.99) / 1.0e3, emax / 1.0e3, γmax, #hide
    ) #hide
end #hide
Markdown.parse(String(take!(io))) #hide

# ## Electron Phase-Space Density Spectra
#
# In relativistic kinematics, the omnidirectional phase-space density $f(E)$ is related to
# the differential energy distribution $dN/dE$ by
#
# ```math
# f(E) \propto \frac{1}{p\,E_{\mathrm{tot}}} \frac{dN}{dE}
#   \propto \frac{1}{\gamma\sqrt{\gamma^2 - 1}} \frac{dN}{dE},
# ```
#
# which reduces to $f(E) \propto E^{-1/2}\,dN/dE$ in the non-relativistic limit.

function energy_psd(Es_eV; edges_eV)
    nbins = length(edges_eV) - 1
    centers_keV = Vector{Float64}(undef, nbins)
    psd = Vector{Float64}(undef, nbins)
    mc2_eV = Electron.m * c^2 / eV
    N = length(Es_eV)
    for j in 1:nbins
        e_lo, e_hi = edges_eV[j], edges_eV[j + 1]
        ec = sqrt(e_lo * e_hi)
        centers_keV[j] = ec / 1.0e3
        cnt = count(e -> e_lo <= e < e_hi, Es_eV)
        if cnt > 0
            γ = 1.0 + ec / mc2_eV
            phase_vol = γ * sqrt(max(γ^2 - 1.0, 1.0e-30))
            psd[j] = (cnt / (N * (e_hi - e_lo))) / phase_vol
        else
            psd[j] = NaN
        end
    end
    return centers_keV, psd
end

edges_eV = 10.0 .^ range(1.3, 5.2; length = 38)
c_sw, psd_sw = energy_psd(E_sw; edges_eV)
norm_fac = maximum(filter(!isnan, psd_sw))

f3 = Figure(size = (1120, 500), fontsize = 14)

ax_spec = Axis(
    f3[1, 1];
    title = "Compound Electron Acceleration Spectra (N = $N_particles)",
    xlabel = L"E\ [\mathrm{keV}]",
    ylabel = L"\mathrm{Normalized\ Phase\ Space\ Density}\ f(E)",
    xscale = log10, yscale = log10,
    limits = ((0.02, 150.0), (1.0e-10, 5.0)),
)
lines!(
    ax_spec, c_sw, psd_sw ./ norm_fac;
    color = :gray55, linewidth = 2.5, label = "Initial solar wind (50 eV)"
)
palette = [:magenta3, :crimson, :forestgreen, :royalblue, :black]
for ((label, Es), col) in zip(cases, palette)
    cs, psd = energy_psd(Es; edges_eV)
    lines!(ax_spec, cs, psd ./ norm_fac; color = col, linewidth = 2.2, label = label)
end
E_ref = 10.0 .^ range(0.0, 1.8; length = 50)
lines!(
    ax_spec, E_ref, 5.0e-3 .* E_ref .^ (-3.5);
    color = :gray35, linestyle = :dash, linewidth = 1.6, label = L"\propto E^{-3.5}"
)
axislegend(ax_spec; position = :lb, framevisible = false)

ax_nt = Axis(
    f3[1, 2];
    title = "Role of multiple foreshock transients",
    xlabel = L"E\ [\mathrm{keV}]",
    ylabel = L"\mathrm{Normalized\ Phase\ Space\ Density}\ f(E)",
    xscale = log10, yscale = log10,
    limits = ((0.02, 150.0), (1.0e-10, 5.0)),
)
lines!(
    ax_nt, c_sw, psd_sw ./ norm_fac;
    color = :gray55, linewidth = 2.5, label = "Initial solar wind (50 eV)"
)
for ((label, Es), col) in zip(nt_cases, (:darkorange, :royalblue, :black))
    cs, psd = energy_psd(Es; edges_eV)
    lines!(ax_nt, cs, psd ./ norm_fac; color = col, linewidth = 2.2, label = label)
end
axislegend(ax_nt; position = :lb, framevisible = false)

f3 = DisplayAs.PNG(f3) #hide

# Both panels mirror the findings of Shi et al. (2025):
#
# - **Ablation (left)**: SDA alone stalls below $1\,\mathrm{keV}$; adding EST, WW, and MSW
#   boundary scattering progressively extends the tail to $\sim 18\,\mathrm{keV}$; adding
#   Fermi reflection off the approaching transients and $3\times$ core betatron compression
#   breaks the spectrum out past $100\,\mathrm{keV}$ ($\gamma > 1.2$).
# - **Multiple transients (right)**: Increasing $n_t$ from $1 \to 2 \to 3$ steadily extends
#   the high-energy cutoff ($\sim 33 \to 71 \to 111\,\mathrm{keV}$) because electrons that
#   survive or escape one transient are re-trapped by the next.
