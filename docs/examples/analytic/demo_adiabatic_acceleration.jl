# # Adiabatic Acceleration
#
# This example decomposes the energization of a guiding center (GC) into the
# three adiabatic acceleration channels:
#
# 1. **Betatron acceleration**: conservation of ``\mu`` while the particle sees
#    ``|\mathbf{B}|`` change. It has a *local* part, ``\mu\,\partial B/\partial t``,
#    and a *convective* part, ``\mu\,\mathbf{v}_E\cdot\nabla B``, which is the work
#    done by the grad-B drift. The two are one mechanism, not two.
# 2. **Fermi acceleration**: the curvature drift crossing the electric field,
#    ``(m v_\parallel^2/B)\,(\hat{\mathbf{b}}\times\boldsymbol{\kappa})\cdot\mathbf{E}``.
# 3. **Parallel acceleration**: direct work of the parallel electric field,
#    ``q v_\parallel\,(\mathbf{E}\cdot\hat{\mathbf{b}})``.
#
# The fields are generated from a vector potential, so they satisfy
# ``\nabla\cdot\mathbf{B} = 0`` and ``\nabla\times\mathbf{E} =
# -\partial\mathbf{B}/\partial t`` exactly, unlike a pair of independently
# prescribed ``\mathbf{E}`` and ``\mathbf{B}``.
#
# References:
# - [Key Notes: Particle Acceleration](https://henry2004y.github.io/KeyNotes/contents/single.html#particle-acceleration)
# - Northrop, T. G. (1963). *The Adiabatic Motion of Charged Particles*.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
import TestParticle as TP
using LinearAlgebra: norm, normalize, tr
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Physical Setup
#
# Take an axisymmetric mirror described by the vector potential
#
# ```math
# \mathbf{A}(\mathbf{x}, t) = \frac{B_0(t)\,g(z)}{2}
#   \begin{bmatrix} -y \\ x \\ 0 \end{bmatrix},
#   \quad g(z) = 1 + \frac{z^2}{L^2}, \quad B_0(t) = B_{00}(1 + \beta t),
# ```
#
# which gives ``\mathbf{B} = \nabla\times\mathbf{A}`` and
#
# ```math
# \mathbf{E} = -\frac{\partial \mathbf{A}}{\partial t} - \nabla\Phi,
#   \quad \Phi = -E_{\parallel 0}\Bigl(z + \frac{z^2}{2 L_E}\Bigr).
# ```
#
# The ``-\partial\mathbf{A}/\partial t`` piece is the azimuthal field induced by the
# compression; ``-\nabla\Phi`` adds a controllable ``E_\parallel`` without breaking
# Faraday's law, since the curl of a gradient vanishes.

## --- Parameters (SI units) ---
const B₀₀ = 1.0e-3   # [T] field strength at the mirror midplane
const L = 30.0       # [m] mirror scale length
const α = 1 / L^2    # [m⁻²]
const β = 0.5        # [1/s] compression rate
const Epar0 = 5.0e-4 # [V/m] parallel electric field at the midplane
const L_E = 30.0     # [m] scale length of the parallel field
const v_mag = 5.0e4  # [m/s] initial speed
const pitch0 = 45.0  # [deg] initial pitch angle

const eV = TP.eV

g(z) = 1 + α * z^2
gp(z) = 2α * z

"""
    make_fields(βv)

The mirror field and its Faraday-consistent electric field for a compression rate
`βv`, using the module-level constants `B₀₀`, `α`, `Epar0` and `L_E`.
"""
function make_fields(βv)
    B(x, t = 0.0) = SA[
        -0.5 * B₀₀ * (1 + βv * t) * gp(x[3]) * x[1],
        -0.5 * B₀₀ * (1 + βv * t) * gp(x[3]) * x[2],
        B₀₀ * (1 + βv * t) * g(x[3]),
    ]
    E(x, t = 0.0) = SA[
        0.5 * B₀₀ * βv * g(x[3]) * x[2],
        -0.5 * B₀₀ * βv * g(x[3]) * x[1],
        Epar0 * (1 + x[3] / L_E),
    ]
    return B, E
end

const mirror_B, mirror_E = make_fields(β)

## --- Species and initial conditions ---
const species = Proton
const x0 = [2L, 0.0, 0.0]
const tspan = (0.0, 0.2)

"Initial full-orbit state with speed `v_mag` at pitch angle `pitch` [deg]."
function init_state(x0v, pitch, B)
    b̂ = normalize(B(x0v, 0.0))
    e1, _ = TP.get_perp_vector(b̂)
    θ = deg2rad(pitch)
    return [x0v..., (v_mag * (cos(θ) * b̂ + sin(θ) * e1))...]
end

## At the midplane the field is along +z, so the pitch angle maps directly onto v.
const stateinit = init_state(x0, pitch0, mirror_B)

# ## Maxwell Check
#
# Before tracing anything, verify that the pair really solves Maxwell's equations.

let x = SA[60.0, 10.0, 15.0], t = 0.05
    Bf, Ef = TP.Field(mirror_B), TP.Field(mirror_E)
    JB = TP.jacobian(Bf, x, t)
    JE = TP.jacobian(Ef, x, t)
    curlE = SA[JE[3, 2] - JE[2, 3], JE[1, 3] - JE[3, 1], JE[2, 1] - JE[1, 2]]
    dBdt = TP.derivative_t(Bf, x, t)

    println("∇⋅B / ‖∇B‖          = ", round(abs(tr(JB)) / norm(JB), sigdigits = 2))
    println("‖∇×E + ∂B/∂t‖/‖∂B/∂t‖ = ", round(norm(curlE + dBdt) / norm(dBdt), sigdigits = 2))
end

# ## GC Tracing
#
# The guiding center is traced with a tight tolerance: the work-rate decomposition
# below is exact along the GC trajectory, so whatever mismatch is left comes from
# the time integration and must be resolved for the check to be meaningful.

stateinit_gc, param_gc = prepare_gc(stateinit, mirror_E, mirror_B; species)
q, q2m, μ, Efunc, Bfunc = param_gc
m = q / q2m

sol_gc = solve(
    ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc), Vern9();
    saveat = range(tspan..., length = 8001), reltol = 1.0e-12, abstol = 1.0e-10
)

vpar = [u[4] for u in sol_gc.u]
nbounce = count(i -> vpar[i] * vpar[i - 1] < 0, 2:length(vpar))
ε = [
    TP.adiabaticity_components(SA[u[1], u[2], u[3]], Bfunc, μ, t; species)
        for (u, t) in zip(sol_gc.u, sol_gc.t)
]
εmax = ntuple(i -> maximum(e[i] for e in ε), 3)

println("Magnetic moment μ = ", round(μ, sigdigits = 4), " J/T")
println("Bounces = ", nbounce)
println("Adiabaticity max over the orbit = ", round.(εmax, sigdigits = 2))

# ## Energy Decomposition
#
# `get_work_rates_gc` returns the four leading-order terms
# ``(P_\parallel, P_\mathrm{Fermi}, P_{\nabla B}, P_{\partial B/\partial t})``.
# Since ``\mathbf{v}_E = (\mathbf{E}\times\hat{\mathbf{b}})/B``,
#
# ```math
# P_{\nabla B} = \frac{\mu}{B}(\hat{\mathbf{b}}\times\nabla B)\cdot\mathbf{E}
#   = \mu\,\mathbf{v}_E\cdot\nabla B,
# ```
#
# so the grad-B drift work is the convective part of the betatron channel rather
# than a fourth mechanism. `get_work_channels` performs that grouping and returns
# ``(P_\parallel, P_\mathrm{Fermi}, P_\mathrm{betatron})``.

"""
    accumulate_work(sol, p)

Instantaneous work rates, their trapezoidal time integrals, the GC kinetic energy
`K = ½ m v∥² + μ B`, and `|B|` along `sol`.
"""
function accumulate_work(sol, p)
    q, q2m, μ, _, Bf = p
    mass = q / q2m
    ts = sol.t
    n = length(ts)
    P = zeros(n, 4)
    K = zeros(n)
    Bmag = zeros(n)

    for (i, t) in enumerate(ts)
        u = sol.u[i]
        X = SA[u[1], u[2], u[3]]
        Bmag[i] = norm(Bf(X, t))
        K[i] = 0.5 * mass * u[4]^2 + μ * Bmag[i]
        P[i, :] .= TP.get_work_rates_gc(u, p, t)
    end

    W = similar(P)
    W[1, :] .= 0
    for i in 2:n
        W[i, :] .= W[i - 1, :] .+ 0.5 .* (P[i, :] .+ P[i - 1, :]) .* (ts[i] - ts[i - 1])
    end

    return ts, P, W, K, Bmag
end

"""
    trace_case(x0v, pitch, tspan; βv = β, nout)

Trace a guiding center launched at `x0v` with pitch angle `pitch` [deg] in a mirror
compressed at rate `βv`, and return `(ts, W, ΔK)` in SI units.
"""
function trace_case(x0v, pitch, tspanv; βv = β, nout = 4001)
    B, E = make_fields(βv)
    sigc, p = prepare_gc(init_state(x0v, pitch, B), E, B; species)
    sol = solve(
        ODEProblem(trace_gc!, sigc, tspanv, p), Vern9();
        saveat = range(tspanv..., length = nout), reltol = 1.0e-12, abstol = 1.0e-10
    )
    ts, _, W, K, _ = accumulate_work(sol, p)

    return ts, W, K .- K[1]
end

ts, P, W, K, Bmag = accumulate_work(sol_gc, param_gc)

## Columns of P and W: 1 = parallel, 2 = Fermi, 3 = grad-B (convective part of
## betatron), 4 = ∂B/∂t (local part of betatron).
W_par = W[:, 1]
W_fermi = W[:, 2]
W_conv = W[:, 3]
W_loc = W[:, 4]
W_beta = W_conv .+ W_loc
W_total = W_par .+ W_fermi .+ W_beta
ΔK = K .- K[1]

using Markdown, Printf #hide
io = IOBuffer() #hide
println(io, "| Quantity | Energy gain [eV] |") #hide
println(io, "| -------- | ---------------- |") #hide
@printf(io, "| ΔK (actual) | %.3f |\n", ΔK[end] / eV) #hide
@printf(io, "| Σ work rates | %.3f |\n", W_total[end] / eV) #hide
@printf(io, "| Betatron | %.3f |\n", W_beta[end] / eV) #hide
@printf(io, "| ├ local μ ∂B/∂t | %.3f |\n", W_loc[end] / eV) #hide
@printf(io, "| └ convective μ v_E⋅∇B | %.3f |\n", W_conv[end] / eV) #hide
@printf(io, "| Fermi | %.3f |\n", W_fermi[end] / eV) #hide
@printf(io, "| E∥ | %.3f |\n", W_par[end] / eV) #hide
Markdown.parse(String(take!(io))) #hide

# ## Visualization

## --- Figure 1: orbit, bounce motion and energy ---
f1 = Figure(size = (1100, 650), fontsize = 14)

ax_orbit = Axis(
    f1[1, 1], title = "GC orbit in the meridional plane",
    xlabel = "ρ [m]", ylabel = "z [m]", aspect = DataAspect(),
)
lines!(
    ax_orbit, [hypot(u[1], u[2]) for u in sol_gc.u], [u[3] for u in sol_gc.u],
    color = :steelblue, linewidth = 1
)

ax_vpar = Axis(
    f1[1, 2], title = "Parallel velocity", xlabel = "t [s]", ylabel = "v∥ [m/s]"
)
lines!(ax_vpar, ts, vpar, color = :coral, linewidth = 1)

ax_B = Axis(
    f1[2, 1], title = "|B| along the orbit", xlabel = "t [s]", ylabel = "|B| [T]"
)
lines!(ax_B, ts, Bmag, color = :steelblue, linewidth = 1)

ax_K = Axis(
    f1[2, 2], title = "K = ½ m v∥² + μB", xlabel = "t [s]", ylabel = "K [eV]"
)
lines!(ax_K, ts, K ./ eV, color = :darkgreen, linewidth = 2)

f1 = DisplayAs.PNG(f1) #hide

# The orbit bounces between the two mirror throats while drifting slowly in radius,
# and ``v_\parallel`` reverses sign once per bounce.

## --- Figure 2: cumulative work per channel ---
f2 = Figure(size = (1000, 800), fontsize = 14)

ax_cum = Axis(
    f2[1, 1], title = "Cumulative work per channel",
    xlabel = "t [s]", ylabel = "Energy gain [eV]",
)
lines!(ax_cum, ts, ΔK ./ eV, color = :black, linewidth = 3, label = L"\Delta K")
lines!(
    ax_cum, ts, W_total ./ eV, color = :black, linestyle = :dash, linewidth = 2,
    label = L"\Sigma\ \mathrm{work\ rates}"
)
lines!(
    ax_cum, ts, W_beta ./ eV, color = :coral, linewidth = 2,
    label = L"\mathrm{Betatron}\ (\mathrm{local}+\mathrm{convective})"
)
lines!(
    ax_cum, ts, W_loc ./ eV, color = :coral, linestyle = :dot, linewidth = 1.5,
    label = L"\mathrm{local}\ \mu\,\partial B/\partial t"
)
lines!(
    ax_cum, ts, W_fermi ./ eV, color = :teal, linewidth = 2, label = L"\mathrm{Fermi}"
)
lines!(
    ax_cum, ts, W_par ./ eV, color = :purple, linewidth = 2, label = L"E_\parallel"
)
axislegend(ax_cum; position = :lt)

ax_inst = Axis(
    f2[2, 1], title = "Instantaneous work rates",
    xlabel = "t [s]", ylabel = "Power [eV/s]",
)
@views lines!(
    ax_inst, ts, P[:, 1] ./ eV, color = :purple, linewidth = 1, label = L"P_\parallel"
)
@views lines!(
    ax_inst, ts, P[:, 2] ./ eV, color = :teal, linewidth = 1, label = L"P_\mathrm{Fermi}"
)
@views lines!(
    ax_inst, ts, P[:, 3] ./ eV, color = :goldenrod, linewidth = 1,
    label = L"P_{\nabla B}\ (\mathrm{convective})"
)
@views lines!(
    ax_inst, ts, P[:, 4] ./ eV, color = :coral, linewidth = 1,
    label = L"P_{\partial B/\partial t}\ (\mathrm{local})"
)
axislegend(ax_inst; position = :rt)

f2 = DisplayAs.PNG(f2) #hide

# Betatron and Fermi both grow secularly and end up within ~20 % of each other.
# The convective part of betatron is *negative*: the azimuthal induced field drives
# an inward ``\mathbf{E}\times\mathbf{B}`` drift, and ``|\mathbf{B}|`` weakens
# toward the axis, so convection cools the perpendicular motion and cancels about
# half of the local ``\mu\,\partial B/\partial t`` gain.

# ## Why the parallel channel stays bounded
#
# The purple curve oscillates and returns to zero instead of accumulating. Part of
# this is structural: the mirror is axisymmetric and ``\mathbf{A}`` is purely
# azimuthal, so the induced field ``-\partial\mathbf{A}/\partial t`` is azimuthal as
# well and therefore perpendicular to the meridional ``\hat{\mathbf{b}}``. It carries
# no ``E_\parallel`` at all, and the whole parallel channel comes from
# ``-\nabla\Phi``. For the parameters used here the bounce motion together with the
# slow azimuthal drift keeps the guiding center close to where it started, so this
# conservative contribution is predominantly oscillatory: the energy taken from
# ``\Phi`` on the way up a bounce is handed back on the way down, and the net
# ``E_\parallel`` work stays near zero.
#
# That is a statement about these parameters, not a theorem. The parallel channel
# integrates ``q\,E_\parallel\,v_\parallel`` along the path the guiding center
# actually samples, and the azimuthal ``\mathbf{E}\times\mathbf{B}`` drift does move
# the particle across equipotential surfaces of ``\Phi``. A steeper or less symmetric
# potential, or a drift that moves systematically in one direction, does let a
# trapped particle gain energy from ``E_\parallel``.
#
# Net parallel acceleration is easiest to see on an orbit that does not return. Launching the same
# particle at the throat with a small pitch angle makes it stream through the mirror
# instead of bouncing, and the channel switches on:

ts_tr, W_tr, ΔK_tr = trace_case([2L, 0.0, -60.0], 10.0, (0.0, 0.004); nout = 4001)

f3 = Figure(size = (900, 500), fontsize = 14)
ax_tr = Axis(
    f3[1, 1], title = "Transiting orbit: cumulative work per channel",
    xlabel = "t [s]", ylabel = "Energy gain [eV]",
)
lines!(ax_tr, ts_tr, ΔK_tr ./ eV, color = :black, linewidth = 3, label = L"\Delta K")
lines!(
    ax_tr, ts_tr, (W_tr[:, 3] .+ W_tr[:, 4]) ./ eV, color = :coral, linewidth = 2,
    label = L"\mathrm{Betatron}"
)
lines!(
    ax_tr, ts_tr, W_tr[:, 2] ./ eV, color = :teal, linewidth = 2, label = L"\mathrm{Fermi}"
)
lines!(
    ax_tr, ts_tr, W_tr[:, 1] ./ eV, color = :purple, linewidth = 2, label = L"E_\parallel"
)
axislegend(ax_tr; position = :lt)
f3 = DisplayAs.PNG(f3) #hide

# Over a single transit the compression has barely acted (``\beta t = 0.002``), so
# betatron is negligible and the parallel channel dominates: the particle slides
# monotonically down ``\Phi`` and converts ``q\,\Delta\Phi`` into kinetic energy once.
# Reversing the motion would hand that energy back, which is approximately what the
# trapped orbit does over each of its bounces.

# ## Validation
#
# The mirror is axisymmetric with no azimuthal field component, so
# ``\hat{\mathbf{b}}`` stays in the meridional plane while
# ``\nabla\times\hat{\mathbf{b}}`` is purely azimuthal and
# ``\hat{\mathbf{b}}\cdot(\nabla\times\hat{\mathbf{b}}) = 0``, hence
# ``B^*_\parallel = |\mathbf{B}|`` and the four work rates sum to ``dK/dt``
# *exactly*, with no ``\mathcal{O}(\epsilon)`` remainder. The residual below must
# therefore fall off with the output resolution, and it does, as ``h^2``.

let
    println("residual |ΔK - ΣW| / ΔK versus number of output points:")
    for nout in (2001, 8001, 32001)
        sol = solve(
            ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc), Vern9();
            saveat = range(tspan..., length = nout), reltol = 1.0e-12, abstol = 1.0e-10
        )
        _, _, Wn, Kn, _ = accumulate_work(sol, param_gc)
        res = abs((Kn[end] - Kn[1]) - sum(Wn[end, :])) / abs(Kn[end] - Kn[1])
        println("  nout = ", nout, ":  ", round(res, sigdigits = 3))
    end
end

# ## Parameter Study
#
# Raising the compression rate ``\beta`` strengthens the induced electric field and
# the field growth together, so betatron and Fermi scale linearly, while the
# parallel channel, being bounded, does not.

β_cases = [0.0, 0.25, 0.5, 0.75, 1.0]
channels = zeros(length(β_cases), 3)

for (k, βv) in enumerate(β_cases)
    _, Wv, _ = trace_case(x0, pitch0, tspan; βv, nout = 4001)
    channels[k, 1] = Wv[end, 1]                  # parallel
    channels[k, 2] = Wv[end, 2]                  # Fermi
    channels[k, 3] = Wv[end, 3] + Wv[end, 4]     # betatron
end

f4 = Figure(size = (850, 500), fontsize = 14)
ax_scan = Axis(
    f4[1, 1], title = "Channel contributions vs. compression rate",
    xlabel = "β [1/s]", ylabel = "Energy gain [eV]",
)
scatterlines!(
    ax_scan, β_cases, channels[:, 1] ./ eV, color = :purple, linewidth = 2,
    markersize = 12, label = L"E_\parallel"
)
scatterlines!(
    ax_scan, β_cases, channels[:, 2] ./ eV, color = :teal, linewidth = 2,
    markersize = 12, label = L"\mathrm{Fermi}"
)
scatterlines!(
    ax_scan, β_cases, channels[:, 3] ./ eV, color = :coral, linewidth = 2,
    markersize = 12, label = L"\mathrm{Betatron}"
)
axislegend(ax_scan; position = :lt)
f4 = DisplayAs.PNG(f4) #hide

# At ``\beta = 0`` the field is static, no electric field is induced, and both
# betatron and Fermi vanish identically: all that remains is the bounded parallel
# channel. Switching the compression on adds the two channels that can actually
# energize a trapped particle.
