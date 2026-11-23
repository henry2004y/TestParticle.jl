# # [Cosmic Ray Propagation](@id Cosmic-Ray)
#
# This example shows how to trace cosmic rays in 3D magnetohydrodynamic (MHD) turbulence in
# two complementary regimes:
#
# 1. **Magnetostatic limit (``\mathbf{E} = \mathbf{0}``)**: high-energy particles whose speed
#    far exceeds the Alfvén speed (``v \gg V_A``) experience negligible motional electric
#    fields over a scattering time, conserving energy and following
#    [Cohet & Marcowith (2016)](https://doi.org/10.1051/0004-6361/201527376).
# 2. **Relativistic electromagnetic turbulence (``\mathbf{E} \neq \mathbf{0}``)**: retaining
#    the convective electric field ``\mathbf{E} = -\mathbf{V} \times \mathbf{B}`` and solving
#    the relativistic Lorentz equation with `TraceProblem{RelativisticOrbit}` captures
#    stochastic second-order Fermi and betatron acceleration.
#
# For general SI↔dimensionless conversions see [Dimensionless Units and Normalization](@ref).
# Dedicated solver accuracy, gyrophase-conservation, and relativistic pusher comparisons
# (`Vern9`, `Boris`, `MultistepBoris4`, `Vay`, and `HigueraCary`) are covered in
# [Boris Method](@ref), [Solver Accuracy Analysis](@ref), and
# [Relativistic vs Non-relativistic Tracing](@ref).

import DisplayAs #hide
using TestParticle, BorisPushers, StaticArrays, LinearAlgebra, Random, Statistics
using TestParticle: get_gc, get_perp_vector, get_mean_magnitude, get_energy,
    get_work_channels, get_adiabaticity, WrapExtrap
using CairoMakie
CairoMakie.activate!(; type = "png") #hide

# ## 1. Cosmic-ray normalization
#
# A particle of charge ``q``, mass ``m``, velocity ``\mathbf{v}``, Lorentz factor
# ``\gamma = (1 - |\mathbf{v}|^2 / c^2)^{-1/2}``, and proper velocity (momentum per unit
# mass) ``\mathbf{u} = \gamma \mathbf{v}`` propagates in an MHD snapshot with mean field
# ``B_0``, total magnetic field ``\mathbf{B} = \mathbf{B}_0 + \delta\mathbf{B}``, fluid
# velocity ``\mathbf{V} = V_A \mathbf{V}^\prime``, and ideal convective electric field
# ``\mathbf{E} = -\mathbf{V} \times \mathbf{B}``:
# ```math
# \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} = \frac{q}{m} \left( \mathbf{E} + \mathbf{v} \times \mathbf{B} \right),
# \qquad
# \frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = \frac{\mathbf{u}}{\gamma}.
# ```
#
# Normalizing velocity by ``c`` (``\mathbf{u}^\prime = \mathbf{u}/c``), magnetic field by
# ``B_0`` (``\mathbf{B}^\prime = \mathbf{B}/B_0``), electric field by ``c B_0``
# (``\mathbf{E}^\prime = -\beta_A (\mathbf{V}^\prime \times \mathbf{B}^\prime)`` with
# ``\beta_A = V_A / c``), time by ``\Omega_0^{-1} = m / (q B_0)``, and length by
# ``d_0 = c / \Omega_0`` sets ``q = m = c = B_0 = 1`` and gives
# ```math
# \frac{\mathrm{d}\mathbf{u}^\prime}{\mathrm{d}t^\prime} = \mathbf{E}^\prime + \frac{\mathbf{u}^\prime}{\gamma} \times \mathbf{B}^\prime,
# \qquad
# \frac{\mathrm{d}\mathbf{x}^\prime}{\mathrm{d}t^\prime} = \frac{\mathbf{u}^\prime}{\gamma},
# \qquad
# \gamma = \sqrt{1 + |\mathbf{u}^\prime|^2}.
# ```
# In these dimensionless units the local gyroradius and relativistic gyroperiod are
# ```math
# r_L^\prime = \frac{|\mathbf{u}_\perp^\prime|}{|\mathbf{B}^\prime|},
# \qquad
# T_g^\prime = \frac{2\pi \gamma}{|\mathbf{B}^\prime|}.
# ```
#
# When ``\mathbf{E}^\prime = \mathbf{0}``, energy and ``\gamma = \gamma_0`` are conserved.
# Rescaling time by ``\gamma_0`` (`t'' = t' / γ0`) reduces the system to the magnetostatic
# `FullOrbit` equations ``\mathrm{d}\mathbf{v}^\prime/\mathrm{d}t^{\prime\prime} = \mathbf{v}^\prime \times \mathbf{B}^\prime``
# and ``\mathrm{d}\mathbf{x}^\prime/\mathrm{d}t^{\prime\prime} = \mathbf{v}^\prime`` with
# ``q = m = B_0 = 1``, where one unit of distance equals the reference gyroradius at
# ``|\mathbf{v}_\perp^\prime| = 1`` and one gyroperiod takes ``2\pi``. The ratio
# ``r_{L0} / L`` sets how many grid cells span one gyration and therefore how strongly
# magnetic fluctuations scatter the particle.

# ## 2. Synthetic periodic MHD turbulence
#
# We construct a periodic, divergence-free 3D MHD snapshot ``(\mathbf{B}, \mathbf{V})`` from
# discrete Fourier modes ``\mathbf{k} = (2\pi / L)\mathbf{n}`` polarized perpendicular to
# ``\mathbf{k}`` (so ``\nabla \cdot \delta\mathbf{B} = 0`` and ``\nabla \cdot \mathbf{V} = 0``)
# with a Kolmogorov-like amplitude weighting ``|\mathbf{n}|^{-5/6}``. Using an inclusive
# periodic grid ``[0, L]`` and copying the periodic boundary values lets `WrapExtrap()` wrap
# smoothly across all box faces.
#
# Crucially, we normalize ``\mathbf{B}`` and ``\mathbf{V}`` *before* calling
# `build_interpolator` and evaluate the convective electric field
# ``\mathbf{E}(\mathbf{x}) = -\beta_A (\mathbf{V}(\mathbf{x}) \times \mathbf{B}(\mathbf{x}))``
# pointwise *after* interpolation. Computing the cross product pointwise guarantees
# ``\mathbf{E} \cdot \mathbf{B} = 0`` to machine precision everywhere inside each cell,
# whereas interpolating precomputed ``\mathbf{E}`` grid values independently of ``\mathbf{B}``
# would introduce spurious parallel electric fields ``E_\parallel``.

function _apply_periodic_endpoints!(A)
    A[:, end, :, :] .= @view A[:, 1, :, :]
    A[:, :, end, :] .= @view A[:, :, 1, :]
    A[:, :, :, end] .= @view A[:, :, :, 1]
    return A
end

function generate_mhd_turbulence(
        gridx, gridy, gridz, L;
        B0 = 1.0, dB_amp = 0.35, dV_amp = 1.0, n_modes = 16, kmax = 3, seed = 42
    )
    nx, ny, nz = length(gridx) - 1, length(gridy) - 1, length(gridz) - 1
    B = zeros(3, nx + 1, ny + 1, nz + 1)
    V = zeros(3, nx + 1, ny + 1, nz + 1)
    B[3, :, :, :] .= B0

    dB_raw = zeros(3, nx + 1, ny + 1, nz + 1)
    dV_raw = zeros(3, nx + 1, ny + 1, nz + 1)
    rng = MersenneTwister(seed)
    k0 = 2π / L

    for _ in 1:n_modes
        nvec = rand(rng, -kmax:kmax, 3)
        while all(iszero, nvec)
            nvec = rand(rng, -kmax:kmax, 3)
        end
        kvec = (k0) * SVector{3, Float64}(nvec[1], nvec[2], nvec[3])
        e1, e2 = get_perp_vector(kvec)
        weight = inv(norm(nvec)^(5 / 6))

        θB = 2π * rand(rng)
        ϕB = 2π * rand(rng)
        aB = randn(rng) * weight
        polB = cos(θB) * e1 + sin(θB) * e2

        θV = 2π * rand(rng)
        ϕV = 2π * rand(rng)
        aV = randn(rng) * weight
        polV = cos(θV) * e1 + sin(θV) * e2

        for kk in 1:nz, j in 1:ny, i in 1:nx
            s = kvec[1] * gridx[i] + kvec[2] * gridy[j] + kvec[3] * gridz[kk]
            cB = aB * cos(s + ϕB)
            cV = aV * cos(s + ϕV)
            dB_raw[1, i, j, kk] += cB * polB[1]
            dB_raw[2, i, j, kk] += cB * polB[2]
            dB_raw[3, i, j, kk] += cB * polB[3]
            dV_raw[1, i, j, kk] += cV * polV[1]
            dV_raw[2, i, j, kk] += cV * polV[2]
            dV_raw[3, i, j, kk] += cV * polV[3]
        end
    end

    dB_rms = get_mean_magnitude(@view dB_raw[:, 1:nx, 1:ny, 1:nz])
    dV_rms = get_mean_magnitude(@view dV_raw[:, 1:nx, 1:ny, 1:nz])
    B .+= (dB_amp * B0 / dB_rms) .* dB_raw
    V .= (dV_amp / dV_rms) .* dV_raw
    _apply_periodic_endpoints!(B)
    _apply_periodic_endpoints!(V)
    return B, V
end

nx = ny = nz = 24
dx = 0.5
L = nx * dx
gridx = range(0.0, L; length = nx + 1)
gridy = range(0.0, L; length = ny + 1)
gridz = range(0.0, L; length = nz + 1)

B0 = 1.0
B_arr, V_arr = generate_mhd_turbulence(gridx, gridy, gridz, L; B0)

Bfunc = build_interpolator(B_arr ./ B0, (gridx, gridy, gridz); bc = WrapExtrap())
Vfunc = build_interpolator(V_arr, (gridx, gridy, gridz); bc = WrapExtrap())

Efunc = let βA = 0.1, Vfunc = Vfunc, Bfunc = Bfunc
    x -> -βA * (Vfunc(x) × Bfunc(x))
end;

# ## 3. Case 1: Magnetostatic turbulence (`E = 0`)
#
# We first neglect the electric field and trace an ensemble with `MultistepBoris4`, which
# preserves kinetic energy to machine precision over hundreds of gyro-orbits while delivering
# 4th-order gyrophase accuracy.

# ### 3.1 Constant-initial-gyroradius injection and ensemble tracing
#
# To compare spatial transport cleanly across energies, we fix the **initial gyroradius**
# ``r_{L0}`` rather than the total speed. Because ``r_L = |\mathbf{v}_\perp| / |\mathbf{B}|``,
# we scale the injected perpendicular speed by the local field magnitude:
# ```math
# |\mathbf{v}_\perp| = r_{L0} |\mathbf{B}|,
# \qquad
# |\mathbf{v}_\parallel| = |\mathbf{v}_\perp| \frac{\mu_0}{\sqrt{1 - \mu_0^2}},
# ```
# where ``\mu_0 = \cos\theta_0`` is the initial pitch-angle cosine.

param_b = prepare(; B = Bfunc, q = 1.0, m = 1.0)

function make_injection(Lx, Ly, Lz; rL = 1.0)
    function prob_func(prob, ctx)
        r = rand(ctx.rng, 5)
        x = Lx * (0.1 + 0.8 * r[1])
        y = Ly * (0.1 + 0.8 * r[2])
        z = Lz * (0.1 + 0.8 * r[3])
        loc = SA[x, y, z]

        Bvec = prob.p.B(loc, 0.0)
        Bmag = norm(Bvec)
        b0 = Bvec ./ Bmag
        e1, e2 = get_perp_vector(b0)

        ϕ = 2π * r[4]
        μ0 = 0.999 * r[5] # pitch-angle cosine in [0, 1)
        sinα = √(1 - μ0^2)
        sinϕ, cosϕ = sincos(ϕ)

        vperp_mag = rL * Bmag
        vperp = (e1 .* cosϕ .+ e2 .* sinϕ) .* vperp_mag
        vpar_mag = vperp_mag * μ0 / sinα
        v0 = vpar_mag .* b0 .+ vperp

        return remake(prob; u0 = SA[x, y, z, v0[1], v0[2], v0[3]])
    end
    return prob_func
end

rL0 = 1.0
stateinit = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
tspan_b = (0.0, 2π * 150)
prob_b = TraceProblem(stateinit, tspan_b, param_b; prob_func = make_injection(L, L, L; rL = rL0))

alg_b = MultistepBoris4(; n = 4)
dt_b = 2π / 40
sols_b = solve(
    prob_b, alg_b, EnsembleThreads();
    dt = dt_b, saveat = 15 * dt_b, trajectories = 32, seed = 1234
);

# ### 3.2 Trajectories, guiding centers, and phase-space evolution
#
# We plot four sample orbits together with their guiding centers (via `get_gc`) and inspect
# the time evolution of the instantaneous gyroradius ``r_L`` and pitch-angle cosine ``\mu``
# for the first trajectory. As the particle traverses varying ``|\mathbf{B}|`` and pitch-angle
# scatters off magnetic fluctuations, ``r_L`` fluctuates around ``r_{L0}`` and ``\mu``
# undergoes random reversals.

f1 = Figure(; fontsize = 22, size = (860, 680))
ax1 = Axis3(
    f1[1, 1];
    title = "Sample trajectories and guiding centers",
    xlabel = L"x",
    ylabel = L"y",
    zlabel = L"z",
    aspect = (1.0, 1.0, 1.25),
    xticks = LinearTicks(4),
    yticks = LinearTicks(4),
    zticks = LinearTicks(5),
    xlabeloffset = 45,
    ylabeloffset = 45,
    zlabeloffset = 65
)

palette = Makie.wong_colors()
for (i, sol) in enumerate(sols_b.u[1:4])
    x = sol[1, :]
    y = sol[2, :]
    z = sol[3, :]
    lines!(ax1, x, y, z; color = palette[i], linewidth = 1.8)
    gcs = [get_gc(u, param_b) for u in sol.u]
    lines!(
        ax1,
        [g[1] for g in gcs],
        [g[2] for g in gcs],
        [g[3] for g in gcs];
        color = palette[i],
        linestyle = :dash,
        linewidth = 1.8
    )
end
Legend(
    f1[0, 1],
    [
        LineElement(; color = :black, linewidth = 2.0),
        LineElement(; color = :black, linestyle = :dash, linewidth = 2.0),
    ],
    ["Full orbit", "Guiding center"];
    orientation = :horizontal,
    framevisible = false,
    tellwidth = false
)
f1 = DisplayAs.PNG(f1) #hide

# Next we inspect how the instantaneous gyroradius ``r_L`` and pitch-angle cosine ``\mu``
# evolve along the first trajectory.

function phase_evolution(sol, Bfunc)
    n = length(sol.u)
    ρ = zeros(n)
    μ = zeros(n)
    for k in 1:n
        x = sol.u[k][SA[1, 2, 3]]
        v = sol.u[k][SA[4, 5, 6]]
        B = Bfunc(x)
        Bmag = norm(B)
        b̂ = B ./ Bmag
        vpar = b̂ ⋅ v
        vperp = v - vpar * b̂
        ρ[k] = norm(vperp) / Bmag
        μ[k] = vpar / norm(v)
    end
    return sol.t, ρ, μ
end

t1, ρ1, μ1 = phase_evolution(sols_b.u[1], Bfunc)

f2 = Figure(; fontsize = 22, size = (1000, 680))
axρ = Axis(f2[1, 1]; ylabel = L"r_L", title = "Single-particle gyroradius and pitch-angle evolution")
axμ = Axis(f2[2, 1]; xlabel = L"t / (2\pi)", ylabel = L"\mu")
linkxaxes!(axρ, axμ)
hidexdecorations!(axρ; grid = false)
ylims!(axρ, 0.0, 1.65)
ylims!(axμ, -1.1, 1.1)

lines!(axρ, t1 ./ (2π), ρ1; color = :steelblue, linewidth = 2.0, label = L"r_L")
hlines!(axρ, [rL0]; color = :tomato, linestyle = :dash, linewidth = 2.0, label = L"r_{L0}")
lines!(axμ, t1 ./ (2π), μ1; color = :seagreen, linewidth = 2.0, label = L"\mu")

axislegend(axρ; position = :rt, framevisible = false)
axislegend(axμ; position = :rt, framevisible = false)
f2 = DisplayAs.PNG(f2) #hide

# ### 3.3 Energy scan: spatial diffusion vs. initial gyroradius
#
# Sweeping ``r_{L0}`` demonstrates that higher-energy particles sample larger turbulent scales
# per gyration and exhibit faster growth in ensemble mean squared displacement
# ``\langle |\Delta\mathbf{x}|^2 \rangle``.

rL_list = [0.5, 1.0, 2.0]

f3 = Figure(; fontsize = 22, size = (1000, 600))
ax3 = Axis(
    f3[1, 1];
    xlabel = L"t / (2\pi)",
    ylabel = L"\langle |\Delta \mathbf{x}|^2 \rangle",
    title = "Mean squared displacement vs. initial gyroradius"
)

for (j, rL) in enumerate(rL_list)
    pfunc = make_injection(L, L, L; rL)
    pscan = TraceProblem(stateinit, (0.0, 2π * 100), param_b; prob_func = pfunc)
    ssols = solve(
        pscan, alg_b, EnsembleThreads();
        dt = dt_b, saveat = 10 * dt_b, trajectories = 32, seed = 1234
    )
    nt = length(ssols.u[1].u)
    msd = zeros(nt)
    for sol in ssols.u
        r0 = sol.u[1][SA[1, 2, 3]]
        for k in 1:nt
            Δ = sol.u[k][SA[1, 2, 3]] - r0
            msd[k] += Δ ⋅ Δ
        end
    end
    msd ./= length(ssols.u)
    lines!(ax3, ssols.u[1].t ./ (2π), msd; color = palette[j], linewidth = 2.5, label = L"r_{L0} = %$rL")
end
axislegend(ax3; position = :lt, framevisible = false)
f3 = DisplayAs.PNG(f3) #hide

# ## 4. Case 2: Relativistic MHD EM turbulence (`E ≠ 0`)
#
# When the convective MHD electric field ``\mathbf{E} = -\beta_A (\mathbf{V} \times \mathbf{B})``
# is retained, particles exchange energy with moving magnetized fluid eddies. We now solve
# the full relativistic Lorentz equation using `TraceProblem{RelativisticOrbit}`, whose
# state stores ``(\mathbf{x}, \mathbf{u})`` with proper velocity ``\mathbf{u} = \gamma\mathbf{v}``.

# ### 4.1 Relativistic setup, mono-energetic injection, and adaptive Higuera–Cary tracing
#
# We prepare both the electromagnetic problem (`E = Efunc`) and a relativistic `E = 0`
# reference on the identical magnetic field (`c = 1`). All particles are injected isotropically
# with a single initial Lorentz factor ``\gamma_0 = 2.0`` (``|\mathbf{u}_0| = \sqrt{\gamma_0^2 - 1}``)
# so any broadening in ``\gamma`` directly measures energization by ``\mathbf{E}``.
#
# Because cosmic-ray gyroperiods scale as ``T_g \propto \gamma / |\mathbf{B}|`` while particles
# gain or lose energy, we integrate with `AdaptiveHigueraCary(; safety = 0.025)`, which is
# volume-preserving, preserves the exact relativistic ``\mathbf{E} \times \mathbf{B}`` drift,
# tracks the local relativistic gyroperiod automatically, and conserves energy to machine
# precision at step nodes when ``\mathbf{E} = \mathbf{0}``.

param_em = prepare(; E = Efunc, B = Bfunc, q = 1.0, m = 1.0, c = 1.0)
param_b_rel = prepare(; B = Bfunc, q = 1.0, m = 1.0, c = 1.0)

function make_relativistic_injection(Lx, Ly, Lz; γ0 = 2.0)
    umag = sqrt(γ0^2 - 1)
    function prob_func(prob, ctx)
        r = rand(ctx.rng, 5)
        x = Lx * (0.1 + 0.8 * r[1])
        y = Ly * (0.1 + 0.8 * r[2])
        z = Lz * (0.1 + 0.8 * r[3])
        loc = SA[x, y, z]

        b0 = normalize(prob.p.B(loc, 0.0))
        e1, e2 = get_perp_vector(b0)
        μ0 = 2 * r[4] - 1
        ϕ = 2π * r[5]
        sinα = √(1 - μ0^2)
        sinϕ, cosϕ = sincos(ϕ)
        uvec = umag * (μ0 * b0 + sinα * (cosϕ * e1 + sinϕ * e2))

        return remake(prob; u0 = SA[x, y, z, uvec[1], uvec[2], uvec[3]])
    end
    return prob_func
end

γ0 = 2.0
u0_mag = sqrt(γ0^2 - 1)
Tg0 = 2π * γ0
tspan_em = (0.0, 100 * Tg0)
saveat_em = range(0.0, 100 * Tg0; length = 201)
stateinit_rel = SA[0.0, 0.0, 0.0, u0_mag, 0.0, 0.0]

pf_rel = make_relativistic_injection(L, L, L; γ0)
prob_em = TraceProblem{RelativisticOrbit}(stateinit_rel, tspan_em, param_em; prob_func = pf_rel)
prob_b_rel = TraceProblem{RelativisticOrbit}(stateinit_rel, tspan_em, param_b_rel; prob_func = pf_rel)

alg_rel = AdaptiveHigueraCary(; safety = 0.025)
sols_em = solve(
    prob_em, alg_rel, EnsembleThreads();
    saveat = saveat_em, trajectories = 64, seed = 1234
)
sols_b_rel = solve(
    prob_b_rel, alg_rel, EnsembleThreads();
    saveat = saveat_em, trajectories = 64, seed = 1234
);

# With `q = m = c = 1`, `get_energy(sol)` returns the kinetic energy ``\gamma - 1``, so
# `1 .+ get_energy(sol)` directly gives ``\gamma(t)``.

γ_em_mat = reduce(hcat, [1.0 .+ get_energy(sol) for sol in sols_em.u])
γ_b_mat = reduce(hcat, [1.0 .+ get_energy(sol) for sol in sols_b_rel.u])
γ_end_em = @view γ_em_mat[end, :]

println("E = 0 final max |γ - γ0| = ", maximum(abs.(γ_b_mat[end, :] .- γ0)))
println(
    "E ≠ 0 final γ: mean = ", round(mean(γ_end_em); digits = 3),
    ", std = ", round(std(γ_end_em); digits = 3)
)

# ### 4.2 Stochastic acceleration and momentum-space diffusion
#
# In the presence of the convective electric field, individual trajectories undergo random
# gains and losses in Lorentz factor ``\gamma(t)`` as they scatter across regions of
# converging or diverging plasma flow. Across the ensemble, the energy variance
# ``\langle (\Delta\gamma)^2 \rangle`` grows roughly diffusively and the initially
# mono-energetic distribution broadens into a wider spectrum with net mean energization,
# while the ``\mathbf{E} = \mathbf{0}`` baseline stays fixed at ``\gamma = \gamma_0``.

t_norm = collect(saveat_em) ./ Tg0
var_γ_em = vec(mean((γ_em_mat .- γ0) .^ 2; dims = 2))
var_γ_b = vec(mean((γ_b_mat .- γ0) .^ 2; dims = 2))

f4 = Figure(; fontsize = 22, size = (1150, 780))

ax_gam = Axis(
    f4[1, 1:2];
    xlabel = L"t / T_{g0}",
    ylabel = L"\gamma",
    title = "Relativistic Lorentz factor evolution"
)
for i in 1:6
    lines!(ax_gam, t_norm, γ_em_mat[:, i]; color = palette[i], linewidth = 2.0)
end
hlines!(
    ax_gam, [γ0];
    color = :black, linestyle = :dash, linewidth = 2.4,
    label = L"\mathbf{E} = \mathbf{0}\ (\gamma = \gamma_0)"
)
axislegend(ax_gam; position = :lt, framevisible = false)

ax_var = Axis(
    f4[2, 1];
    xlabel = L"t / T_{g0}",
    ylabel = L"\langle (\Delta \gamma)^2 \rangle",
    title = "Energy-spread growth"
)
lines!(ax_var, t_norm, var_γ_em; color = :steelblue, linewidth = 2.5, label = L"\mathbf{E} \neq \mathbf{0}")
lines!(ax_var, t_norm, var_γ_b; color = :black, linestyle = :dash, linewidth = 2.2, label = L"\mathbf{E} = \mathbf{0}")
axislegend(ax_var; position = :lt, framevisible = false)

ax_hist = Axis(
    f4[2, 2];
    xlabel = L"\gamma / \gamma_0",
    ylabel = "Count",
    title = "Final energy distribution"
)
hist!(
    ax_hist, γ_end_em ./ γ0;
    bins = 14, color = (:steelblue, 0.75), strokecolor = :white, strokewidth = 1.0
)
vlines!(ax_hist, [1.0]; color = :black, linestyle = :dash, linewidth = 2.4, label = L"\gamma = \gamma_0")
vlines!(
    ax_hist, [mean(γ_end_em) / γ0];
    color = :tomato, linestyle = :dot, linewidth = 3.0, label = L"\langle \gamma \rangle / \gamma_0"
)
axislegend(ax_hist; position = :rt, framevisible = false)

rowgap!(f4.layout, 18)
colgap!(f4.layout, 24)
f4 = DisplayAs.PNG(f4) #hide

# ### 4.3 Acceleration channels and spatial transport
#
# Using `get_work_channels` with `model = RelativisticOrbit()`, we evaluate the leading-order
# adiabatic energization channels at every saved step: parallel electric acceleration
# ``P_\parallel = q v_\parallel (\mathbf{E} \cdot \hat{\mathbf{b}})``, curvature-drift Fermi
# work ``P_{\mathrm{Fermi}} = \gamma m v_\parallel^2 (\hat{\mathbf{b}} \times \boldsymbol{\kappa}) \cdot \mathbf{E} / |\mathbf{B}|``,
# and betatron/grad-``B`` work ``P_{\mathrm{betatron}} = (\mu / |\mathbf{B}|)(\hat{\mathbf{b}} \times \nabla B) \cdot \mathbf{E} + \mu \partial_t B``.
#
# Three physical features stand out:
# 1. Because ``\mathbf{E}(\mathbf{x}) = -\beta_A (\mathbf{V}(\mathbf{x}) \times \mathbf{B}(\mathbf{x}))``
#    is evaluated pointwise after interpolation, ``W_\parallel`` vanishes to machine precision.
# 2. Curvature-drift Fermi heating (``W_{\mathrm{Fermi}} > 0``) dominates over grad-``B``
#    cooling (``W_{\mathrm{betatron}} < 0``).
# 3. In this turbulence, cosmic rays have finite adiabaticity ``\epsilon = r_L / R_c``
#    (computed via `get_adiabaticity`), so the leading-order guiding-center sum
#    ``W_{\mathrm{GC}} = W_\parallel + W_{\mathrm{Fermi}} + W_{\mathrm{betatron}}`` captures
#    only part of the total energy gain ``\langle \Delta\gamma \rangle``; the excess comes
#    from non-adiabatic gyro-resonant scattering, which requires full-orbit
#    `RelativisticOrbit` integration rather than a guiding-center approximation.

nt_em = length(saveat_em)
P_mean = zeros(3, nt_em)
for sol in sols_em.u
    for k in 1:nt_em
        wc = get_work_channels(sol.u[k], param_em, sol.t[k]; model = RelativisticOrbit())
        P_mean[1, k] += wc[1]
        P_mean[2, k] += wc[2]
        P_mean[3, k] += wc[3]
    end
end
P_mean ./= length(sols_em.u)

dt_save = step(saveat_em)
W_mean = zeros(3, nt_em)
for k in 2:nt_em, c in 1:3
    W_mean[c, k] = W_mean[c, k - 1] + 0.5 * (P_mean[c, k - 1] + P_mean[c, k]) * dt_save
end

W_par = @view W_mean[1, :]
W_fermi = @view W_mean[2, :]
W_betatron = @view W_mean[3, :]
W_gc = W_par .+ W_fermi .+ W_betatron
Δγ_mean = vec(mean(γ_em_mat .- γ0; dims = 2))

ϵ0 = mean(
    get_adiabaticity(sol.u[1], param_em, sol.t[1]; model = RelativisticOrbit()) for
        sol in sols_em.u
)
ϵ1 = mean(
    get_adiabaticity(sol.u[end], param_em, sol.t[end]; model = RelativisticOrbit()) for
        sol in sols_em.u
)
println("max |W_par| = ", maximum(abs, W_par))
println("Mean adiabaticity ϵ: initial = ", round(ϵ0; digits = 3), ", final = ", round(ϵ1; digits = 3))

# We now plot the cumulative guiding-center work channels alongside ``\langle \Delta\gamma \rangle``
# and compare spatial diffusion with and without the convective electric field.

function compute_msd(sols, nt)
    msd = zeros(nt)
    for sol in sols.u
        r0 = sol.u[1][SA[1, 2, 3]]
        for k in 1:nt
            Δ = sol.u[k][SA[1, 2, 3]] - r0
            msd[k] += Δ ⋅ Δ
        end
    end
    return msd ./ length(sols.u)
end

msd_em = compute_msd(sols_em, nt_em)
msd_b_rel = compute_msd(sols_b_rel, nt_em)

f5 = Figure(; fontsize = 22, size = (1200, 540))

ax_work = Axis(
    f5[1, 1];
    xlabel = L"t / T_{g0}",
    ylabel = L"\langle \Delta\gamma \rangle,\; \langle W \rangle",
    title = "Cumulative relativistic work channels"
)
lines!(ax_work, t_norm, W_fermi; color = palette[1], linewidth = 2.5, label = L"W_{\mathrm{Fermi}}")
lines!(ax_work, t_norm, W_betatron; color = palette[2], linewidth = 2.5, label = L"W_{\mathrm{betatron}}")
lines!(ax_work, t_norm, W_par; color = palette[3], linewidth = 2.2, label = L"W_\parallel")
lines!(ax_work, t_norm, W_gc; color = palette[4], linewidth = 2.5, label = L"W_{\mathrm{GC}}")
lines!(
    ax_work, t_norm, Δγ_mean;
    color = :black, linestyle = :dash, linewidth = 2.4, label = L"\langle \Delta\gamma \rangle"
)
axislegend(ax_work; position = :lt, nbanks = 2, framevisible = false)

ax_msd = Axis(
    f5[1, 2];
    xlabel = L"t / T_{g0}",
    ylabel = L"\langle |\Delta \mathbf{x}|^2 \rangle",
    title = "Spatial diffusion comparison"
)
lines!(ax_msd, t_norm, msd_em; color = :steelblue, linewidth = 2.5, label = L"\mathbf{E} \neq \mathbf{0}")
lines!(ax_msd, t_norm, msd_b_rel; color = :black, linestyle = :dash, linewidth = 2.4, label = L"\mathbf{E} = \mathbf{0}")
axislegend(ax_msd; position = :lt, framevisible = false)

colgap!(f5.layout, 24)
f5 = DisplayAs.PNG(f5) #hide

# ## 5. Production notes
#
# Moving from this synthetic example to production runs on 3D MHD simulation snapshots is
# straightforward:
#
# - **Field loading and interpolation**: Load the 3D magnetic field `B` and fluid velocity
#   `V` (plus current density `J` if including a resistive term ``\eta \mathbf{J}``),
#   normalize them by ``B_0`` and ``V_A`` on an inclusive periodic grid, and wrap each with
#   `build_interpolator(arr, (gridx, gridy, gridz); bc = WrapExtrap())`.
# - **Electric field evaluation**: Construct ``\mathbf{E}(\mathbf{x}) = -\beta_A (\mathbf{V}(\mathbf{x}) \times \mathbf{B}(\mathbf{x}))``
#   (and optional resistive ``\eta \mathbf{J}(\mathbf{x})``) pointwise after spatial
#   interpolation to avoid interpolation-induced parallel electric fields.
# - **Solver selection**: Use `MultistepBoris4` for magnetostatic (``\mathbf{E} = \mathbf{0}``)
#   propagation, `AdaptiveHigueraCary` on CPUs when ``\gamma`` varies dynamically in
#   electromagnetic turbulence, or fixed-step `HigueraCary` with `EnsembleKernel` for GPU
#   runs (see `demo_gpu`).
# - **Diagnostics and storage**: Compute guiding centers (`get_gc`), energies (`get_energy`),
#   gyroradii (`get_gyroradius`), adiabatic work channels (`get_work_channels`), and field
#   geometry (`get_magnetic_properties`), and persist ensembles to JLD2 via `save_trajectory`.
