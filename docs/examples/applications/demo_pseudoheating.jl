# # Ion Pseudoheating by Low-Frequency Alfvén Waves
#
# This example demonstrates **ion pseudoheating** by low-frequency Alfvén waves, following
# [Dong & Singh (2013), *Ion pseudoheating by low-frequency Alfvén waves revisited*](https://doi.org/10.1063/1.4789608)
# (see also [Dong (2014)](https://doi.org/10.1063/1.4863833) for minor ions).
#
# When low-frequency Alfvén waves ($\omega \ll \Omega_p$) propagate along a background
# magnetic field $\mathbf{B}_0 = B_0 \hat{\mathbf{z}}$, the wave electric field drives a
# transverse $\mathbf{E} \times \mathbf{B}$ drift
#
# ```math
# \mathbf{v}_E(z, t) = \frac{\mathbf{E}_W \times \mathbf{B}_0}{B_0^2}
#   = -\frac{v_A}{B_0}\,\mathbf{B}_W(z, t),
# ```
#
# which oscillates in antiphase with the wave magnetic perturbation $\mathbf{B}_W$. If an
# observation or diagnostic averages particles over a spatial scale larger than the Alfvén
# wavelength ($\Delta z \gg \lambda_A$, low spatial resolution), this spatially varying bulk
# drift broadens the measured velocity distribution function (VDF) and raises the **apparent
# perpendicular temperature** to
#
# ```math
# \frac{T_{p\perp}(t)}{T_0} \simeq 1 + \frac{1}{\beta_p}\,\frac{B_W^2(t)}{B_0^2},
# \qquad \beta_p = \frac{2 v_p^2}{v_A^2},
# ```
#
# where $v_p = \sqrt{k_B T_0 / m_p}$ is the 1D initial thermal speed. Unlike genuine
# thermodynamic heating, pseudoheating is **reversible** and non-dissipative: the magnetic
# moment is conserved, the apparent temperature returns to $T_0$ once the waves subside, and
# a **locally resolved** VDF ($\Delta z \ll \lambda_A$) is simply a shifted copy of the
# initial unheated Maxwellian.

import DisplayAs #hide
using TestParticle, BorisPushers, StaticArrays
using Random, Statistics, Printf
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Dimensionless Alfvén Wave Fields
#
# We work in dimensionless cyclotron–Alfvén units ($q = m = B_0 = v_A = \Omega_p = 1$), so
# time is measured in $\Omega_p^{-1}$, length in $v_A \Omega_p^{-1}$, and velocity in $v_A$.
# A spectrum of $N$ equal-amplitude, left-hand circularly polarized Alfvén waves propagating
# along $\hat{\mathbf{z}}$ has fields
#
# ```math
# \mathbf{B}_W(z, t) = \sqrt{\frac{\epsilon(t)}{N}}\,B_0
#   \sum_{j=1}^{N} \bigl(\cos\phi_j\,\hat{\mathbf{x}} - \sin\phi_j\,\hat{\mathbf{y}}\bigr),
# \qquad
# \mathbf{E}_W(z, t) = -v_A \hat{\mathbf{z}} \times \mathbf{B}_W(z, t),
# ```
#
# where $\phi_j = k_j(v_A t - z) + \varphi_j$, $k_j = \omega_j / v_A$, and the wave energy
# density $\epsilon(t) = B_W^2(t) / B_0^2$ ramps smoothly up to $\epsilon_0 = 0.05$, stays
# flat for $t_1 \le \Omega_p t \le t_2$, and then decays back to zero:
#
# ```math
# \epsilon(t) = \begin{cases}
#   \epsilon_0 e^{-(t - t_1)^2 / \tau^2}, & t < t_1, \\
#   \epsilon_0, & t_1 \le t \le t_2, \\
#   \epsilon_0 e^{-(t - t_2)^2 / \tau^2}, & t > t_2.
# \end{cases}
# ```

const t₁ = 500.0
const t₂ = 1000.0
const τ = 200.0
const ϵ₀ = 0.05
const t_end = 1500.0
const v_p = 0.07
const β_p = 2 * v_p^2
const L_z = 3000.0

@inline function wave_envelope(t)
    if t < t₁
        return ϵ₀ * exp(-((t - t₁) / τ)^2)
    elseif t <= t₂
        return ϵ₀
    else
        return ϵ₀ * exp(-((t - t₂) / τ)^2)
    end
end

struct AlfvenWaveSpectrum{N}
    k₁::Float64
    Δk::Float64
    eiφ::SVector{N, ComplexF64}
    inv_sqrtN::Float64
end

function AlfvenWaveSpectrum(k₁::Float64, k_end::Float64, φs::AbstractVector{Float64})
    N = length(φs)
    Δk = N > 1 ? (k_end - k₁) / (N - 1) : 0.0
    eiφ = SVector{N, ComplexF64}(cis.(φs))
    return AlfvenWaveSpectrum{N}(k₁, Δk, eiφ, inv(sqrt(Float64(N))))
end

@inline function wave_Bperp(w::AlfvenWaveSpectrum{N}, z, t) where {N}
    amp = sqrt(wave_envelope(t)) * w.inv_sqrtN
    ξ = t - z
    z_cur = cis(w.k₁ * ξ)
    z_step = cis(w.Δk * ξ)
    acc = ComplexF64(0.0, 0.0)
    @inbounds for j in 1:N
        acc += z_cur * w.eiφ[j]
        z_cur *= z_step
    end
    return amp * real(acc), -amp * imag(acc)
end

struct AlfvenB{N}
    w::AlfvenWaveSpectrum{N}
end
@inline function (f::AlfvenB)(x, t)
    bx, by = wave_Bperp(f.w, x[3], t)
    return SA[bx, by, 1.0]
end

struct AlfvenE{N}
    w::AlfvenWaveSpectrum{N}
end
@inline function (f::AlfvenE)(x, t)
    bx, by = wave_Bperp(f.w, x[3], t)
    return SA[by, -bx, 0.0]
end

# ## Tracing Proton Ensembles with the Boris Pusher
#
# Following Dong & Singh (2013), we compare two wave configurations with the same total
# energy $\epsilon_0 = 0.05$:
#
# 1. **Monochromatic wave** ($N = 1$, $\omega = 0.05\,\Omega_p$),
# 2. **Broadband wave spectrum** ($N = 41$, $\omega_j \in [0.01, 0.05]\,\Omega_p$, random
#    phases $\varphi_j$).
#
# Protons are initialized uniformly in $z \in [0, 3000]\,v_A \Omega_p^{-1}$ with an
# isotropic Maxwellian velocity distribution ($v_p = 0.07\,v_A$) and integrated with the
# energy-conserving `Boris()` pusher (`dt = 0.25`, ~25 steps per gyroperiod).

function run_pseudoheating(
        w::AlfvenWaveSpectrum; N_particles::Int = 15_000, dt::Float64 = 0.25, seed::Int = 42
    )
    param = prepare(; E = AlfvenE(w), B = AlfvenB(w), q = 1.0, m = 1.0)
    tsave = collect(0.0:25.0:t_end)
    function prob_func(prob, ctx)
        rng = ctx.rng
        z0 = L_z * rand(rng)
        vx0 = v_p * randn(rng)
        vy0 = v_p * randn(rng)
        vz0 = v_p * randn(rng)
        return remake(prob; u0 = SA[0.0, 0.0, z0, vx0, vy0, vz0])
    end
    u0 = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    prob = TraceProblem(u0, (0.0, t_end), param; prob_func)
    sols = solve(
        prob, Boris(), EnsembleThreads();
        dt, saveat = tsave, trajectories = N_particles, seed,
    )
    return sols, tsave
end

w_mono = AlfvenWaveSpectrum(0.05, 0.05, [0.0])
w_spec = AlfvenWaveSpectrum(0.01, 0.05, 2π .* rand(Xoshiro(1234), 41))

sols_mono, tsave = run_pseudoheating(w_mono);
sols_spec, _ = run_pseudoheating(w_spec);

# ## Apparent Perpendicular Temperature and Reversibility
#
# We compute the global apparent perpendicular temperature
# $T_{p\perp}(t) / T_0 = \langle |\mathbf{v}_\perp - \langle \mathbf{v}_\perp \rangle|^2 \rangle / (2 v_p^2)$
# across all particles and compare it with $1 + \epsilon(t)/\beta_p$.

function apparent_temp_ratio(sols, nt)
    N = length(sols.u)
    Tratio = Vector{Float64}(undef, nt)
    for it in 1:nt
        sx = 0.0
        sy = 0.0
        for sol in sols.u
            u = sol.u[it]
            sx += u[4]
            sy += u[5]
        end
        mx = sx / N
        my = sy / N
        v2 = 0.0
        for sol in sols.u
            u = sol.u[it]
            v2 += (u[4] - mx)^2 + (u[5] - my)^2
        end
        Tratio[it] = (v2 / N) / (2 * v_p^2)
    end
    return Tratio
end

T_mono = apparent_temp_ratio(sols_mono, length(tsave))
T_spec = apparent_temp_ratio(sols_spec, length(tsave))
T_analytic = [1.0 + wave_envelope(t) / β_p for t in tsave]
ϵ_curve = [wave_envelope(t) for t in tsave]

idx600 = findfirst(==(600.0), tsave)
@printf(
    "Peak T_p⊥/T₀ at Ω_p t = 600: mono = %.2f, spectrum = %.2f, analytic = %.2f\n",
    T_mono[idx600], T_spec[idx600], T_analytic[idx600]
)
@printf(
    "Final T_p⊥/T₀ at Ω_p t = 1500: mono = %.2f, spectrum = %.2f\n",
    T_mono[end], T_spec[end]
)

f1 = Figure(size = (820, 420), fontsize = 14)
ax1 = Axis(
    f1[1, 1];
    title = "Apparent perpendicular temperature vs. wave energy density",
    xlabel = L"\Omega_p t",
    ylabel = L"T_{p\perp} / T_0",
    limits = ((0, t_end), (0, 7.2)),
)
ax2 = Axis(
    f1[1, 1];
    yaxisposition = :right,
    ylabel = L"B_W^2 / B_0^2",
    limits = ((0, t_end), (0, 0.062)),
    xgridvisible = false, ygridvisible = false,
)
hidexdecorations!(ax2)

l_mono = lines!(ax1, tsave, T_mono; color = :steelblue, linewidth = 2.2)
l_spec = lines!(ax1, tsave, T_spec; color = :crimson, linewidth = 2.2)
l_ana = lines!(
    ax1, tsave, T_analytic; color = :black, linestyle = :dash, linewidth = 2.0
)
l_bw = lines!(
    ax2, tsave, ϵ_curve; color = :gray45, linestyle = :dot, linewidth = 2.2
)
axislegend(
    ax1,
    [l_mono, l_spec, l_ana, l_bw],
    [
        "Monochromatic (N = 1)",
        "Wave spectrum (N = 41)",
        L"Analytic $1 + \epsilon(t)/\beta_p$",
        L"B_W^2 / B_0^2\ \mathrm{(right)}",
    ];
    position = :rt, framevisible = false,
)
f1 = DisplayAs.PNG(f1) #hide

# Both the monochromatic wave and the 41-mode spectrum closely follow the analytic
# $\mathbf{E} \times \mathbf{B}$ prediction (the slight overshoot for the monochromatic
# $\omega = 0.05\,\Omega_p$ wave is the $\mathcal{O}(\omega/\Omega_p)$ polarization-drift
# correction) and return cleanly to $T_{p\perp} / T_0 \approx 1$ as the waves subside.

# ## Phase-Space and Velocity-Space Signatures
#
# At peak wave amplitude ($\Omega_p t = 600$), plotting $(z, v_x)$ phase space alongside
# $-B_{Wx}(z, t) / B_0$ shows that protons track the local $\mathbf{E} \times \mathbf{B}$
# drift velocity $v_{Ex} = -v_A B_{Wx} / B_0$ directly. In $(v_x, v_y)$ velocity space,
# however, the monochromatic wave produces a ring-like distribution whereas the 41-mode
# spectrum fills a smooth disk.

z_mono = [sol.u[idx600][3] for sol in sols_mono.u]
vx_mono = [sol.u[idx600][4] for sol in sols_mono.u]
vy_mono = [sol.u[idx600][5] for sol in sols_mono.u]

z_spec = [sol.u[idx600][3] for sol in sols_spec.u]
vx_spec = [sol.u[idx600][4] for sol in sols_spec.u]
vy_spec = [sol.u[idx600][5] for sol in sols_spec.u]

f2 = Figure(size = (1060, 680), fontsize = 13)

z_grid = range(1000.0, 2400.0; length = 600)
vEx_mono = [-wave_Bperp(w_mono, z, 600.0)[1] for z in z_grid]
vEx_spec = [-wave_Bperp(w_spec, z, 600.0)[1] for z in z_grid]

ax_pm = Axis(
    f2[1, 1];
    title = L"\mathrm{Monochromatic\ wave}\ (N = 1),\ \Omega_p t = 600",
    xlabel = L"z\,\Omega_p / v_A", ylabel = L"v_x / v_A",
    limits = ((1000, 2400), (-0.6, 0.6)),
)
scatter!(
    ax_pm, z_mono, vx_mono;
    color = (:steelblue, 0.25), markersize = 2.5,
)
lines!(
    ax_pm, z_grid, vEx_mono;
    color = :crimson, linewidth = 1.6, label = L"v_{Ex} = -v_A B_{Wx}/B_0",
)
axislegend(ax_pm; position = :rt, framevisible = false)

ax_vm = Axis(
    f2[1, 2];
    title = "Velocity space (N = 1)",
    xlabel = L"v_x / v_A", ylabel = L"v_y / v_A",
    aspect = DataAspect(), limits = ((-0.6, 0.6), (-0.6, 0.6)),
)
scatter!(ax_vm, vx_mono, vy_mono; color = (:steelblue, 0.2), markersize = 2.5)

ax_ps = Axis(
    f2[2, 1];
    title = L"\mathrm{Wave\ spectrum}\ (N = 41),\ \Omega_p t = 600",
    xlabel = L"z\,\Omega_p / v_A", ylabel = L"v_x / v_A",
    limits = ((1000, 2400), (-0.6, 0.6)),
)
scatter!(
    ax_ps, z_spec, vx_spec;
    color = (:steelblue, 0.25), markersize = 2.5,
)
lines!(
    ax_ps, z_grid, vEx_spec;
    color = :crimson, linewidth = 1.6, label = L"v_{Ex} = -v_A B_{Wx}/B_0",
)
vspan!(ax_ps, 1490.0, 1510.0; color = (:darkorange, 0.3))
axislegend(ax_ps; position = :rt, framevisible = false)

ax_vs = Axis(
    f2[2, 2];
    title = "Velocity space (N = 41)",
    xlabel = L"v_x / v_A", ylabel = L"v_y / v_A",
    aspect = DataAspect(), limits = ((-0.6, 0.6), (-0.6, 0.6)),
)
scatter!(ax_vs, vx_spec, vy_spec; color = (:steelblue, 0.2), markersize = 2.5)

colsize!(f2.layout, 1, Relative(0.64))
f2 = DisplayAs.PNG(f2) #hide

# ## Spatial Resolution: Global vs. Local Velocity Distributions
#
# The key observational consequence highlighted by Dong & Singh (2013) is the role of
# **spatial resolution**:
#
# - **Low spatial resolution (global average)**: Under a broadband Alfvén wave spectrum,
#   averaging over many wavelengths produces a broadened, nearly Maxwellian VDF that closely
#   mimics genuine heating (unlike the clearly non-Maxwellian double-peaked VDF of a single
#   monochromatic wave).
# - **High spatial resolution (local window)**: Restricting the sample to a narrow spatial
#   slab ($1490 < z \Omega_p / v_A < 1510$, orange band above) reveals that the local VDF is
#   merely **shifted** by the local $\mathbf{E} \times \mathbf{B}$ drift
#   $\langle v_x \rangle_{\mathrm{local}}$. Transforming into the local mean-velocity frame
#   recovers the original unheated Maxwellian.

function normalized_vdf(vs; edges)
    nbins = length(edges) - 1
    centers = [0.5 * (edges[i] + edges[i + 1]) for i in 1:nbins]
    counts = zeros(Float64, nbins)
    lo, hi = first(edges), last(edges)
    inv_dx = nbins / (hi - lo)
    for v in vs
        if lo <= v < hi
            idx = clamp(floor(Int, (v - lo) * inv_dx) + 1, 1, nbins)
            counts[idx] += 1.0
        end
    end
    m = maximum(counts)
    if m > 0
        counts ./= m
    end
    return centers, counts
end

vx_mono_end = [sol.u[end][4] for sol in sols_mono.u]
vx_spec_end = [sol.u[end][4] for sol in sols_spec.u]

local_mask = 1490.0 .< z_spec .< 1510.0
vx_local = vx_spec[local_mask]
vy_local = vy_spec[local_mask]
vx_local_mean = mean(vx_local)
T_local_ratio = (var(vx_local) + var(vy_local)) / (2 * v_p^2)

@printf(
    "Local slab (1490 < z < 1510, n = %d): ⟨v_x⟩ = %.3f v_A, T_local/T₀ = %.2f\n",
    length(vx_local), vx_local_mean, T_local_ratio
)

edges_glob = range(-0.6, 0.6; length = 61)
edges_loc = range(-0.55, 0.55; length = 27)
v_grid = range(-0.6, 0.6; length = 300)
f_maxw_init = exp.(-(v_grid .^ 2) ./ (2 * v_p^2))

f3 = Figure(size = (1140, 390), fontsize = 13)

ax_a = Axis(
    f3[1, 1];
    title = "(a) Global VDF: Monochromatic (N = 1)",
    xlabel = L"v_x / v_A", ylabel = L"f(v_x) / f_{\max}",
    limits = ((-0.6, 0.6), (0, 1.18)),
)
lines!(
    ax_a, v_grid, f_maxw_init; color = :black, linewidth = 2.0, label = "Initial (t = 0)"
)
c_m600, h_m600 = normalized_vdf(vx_mono; edges = edges_glob)
c_mend, h_mend = normalized_vdf(vx_mono_end; edges = edges_glob)
lines!(
    ax_a, c_m600, h_m600;
    color = :crimson, linestyle = :dash, linewidth = 2.2, label = "Peak (t = 600)",
)
lines!(
    ax_a, c_mend, h_mend;
    color = :steelblue, linestyle = :dot, linewidth = 2.2, label = "Final (t = 1500)",
)
axislegend(ax_a; position = :rt, framevisible = false, labelsize = 11)

ax_b = Axis(
    f3[1, 2];
    title = "(b) Global VDF: Spectrum (N = 41)",
    xlabel = L"v_x / v_A", ylabel = L"f(v_x) / f_{\max}",
    limits = ((-0.6, 0.6), (0, 1.18)),
)
lines!(
    ax_b, v_grid, f_maxw_init; color = :black, linewidth = 2.0, label = "Initial (t = 0)"
)
c_s600, h_s600 = normalized_vdf(vx_spec; edges = edges_glob)
c_send, h_send = normalized_vdf(vx_spec_end; edges = edges_glob)
lines!(
    ax_b, c_s600, h_s600;
    color = :crimson, linestyle = :dash, linewidth = 2.2, label = "Peak (t = 600)",
)
lines!(
    ax_b, c_send, h_send;
    color = :steelblue, linestyle = :dot, linewidth = 2.2, label = "Final (t = 1500)",
)
axislegend(ax_b; position = :rt, framevisible = false, labelsize = 11)

ax_c = Axis(
    f3[1, 3];
    title = "(c) Local VDF (1490 < z < 1510), t = 600",
    xlabel = L"v_x / v_A", ylabel = L"f(v_x) / f_{\max}",
    limits = ((-0.6, 0.6), (0, 1.18)),
)
lines!(
    ax_c, v_grid, f_maxw_init; color = :black, linewidth = 2.0, label = "Initial (t = 0)"
)
c_lab, h_lab = normalized_vdf(vx_local; edges = edges_loc)
c_com, h_com = normalized_vdf(vx_local .- vx_local_mean; edges = edges_loc)
nz_lab = h_lab .> 0
nz_com = h_com .> 0
scatter!(
    ax_c, c_lab[nz_lab], h_lab[nz_lab];
    color = :crimson, markersize = 7, label = "Lab frame",
)
scatter!(
    ax_c, c_com[nz_com], h_com[nz_com];
    color = :white, strokecolor = :steelblue, strokewidth = 1.8, markersize = 7,
    label = "Local mean frame",
)
axislegend(ax_c; position = :lt, framevisible = false, labelsize = 11)

f3 = DisplayAs.PNG(f3) #hide

# In test-particle tracing, the prescribed Alfvén wave fields do not damp, so the apparent
# heating is purely reversible pseudoheating. In self-consistent hybrid or PIC simulations,
# nonlinear wave–particle interactions (such as nonlinear Landau damping and parametric
# instabilities of finite-amplitude Alfvén wave packets) can dissipate part of the wave
# energy into genuine irreversible heating alongside the reversible
# $\mathbf{E} \times \mathbf{B}$ pseudoheating shown here.
