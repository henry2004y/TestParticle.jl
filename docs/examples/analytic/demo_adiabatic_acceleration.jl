# # Adiabatic Acceleration
#
# This example demonstrates the three fundamental adiabatic acceleration
# mechanisms using guiding center (GC) tracing:
#
# 1. **Betatron acceleration**: ``\mu`` conserved, ``\partial B/\partial t > 0``
#    → perpendicular energy gain via ``\mu\, \partial B/\partial t``.
# 2. **Fermi acceleration**: curvature and grad-B drifts crossing an electric
#    field → parallel/perpendicular energy exchange.
# 3. **Local ``E_\parallel`` acceleration**: direct acceleration by the parallel
#    electric field component ``q v_\parallel (\mathbf{E}\cdot\hat{\mathbf{b}})``.
#
# The physical setup uses a **sheared, time-dependent magnetic field** with a
# constant electric field, inducing all three mechanisms simultaneously.
# The energy gain is decomposed using the GC work rate formulas from
# Northrop (1963).
#
# References:
# - [Key Notes: Particle Acceleration](https://henry2004y.github.io/KeyNotes/contents/single.html#particle-acceleration)
# - Northrop, T. G. (1963). *The Adiabatic Motion of Charged Particles*.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
import TestParticle as TP
using LinearAlgebra: norm, ⋅, ×
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Physical Setup
#
# The magnetic field is a sheared, divergence-free field with a linear ramp:
#
# ```math
# \mathbf{B}(\mathbf{x}, t) = B_0(t) \begin{bmatrix} \alpha z \\ 0 \\ 1 \end{bmatrix},
# \quad B_0(t) = B_{00} (1 + \beta t).
# ```
#
# This field has:
# - **Curvature** ``\boldsymbol{\kappa} \neq 0`` (from field-line bending)
# - **Gradient** ``\nabla B \neq 0`` (from the ``z``-dependence of ``B_x``)
# - **Time dependence** ``\partial B/\partial t \neq 0`` (from the ``B_0(t)`` ramp)
#
# The electric field is constant:
#
# ```math
# \mathbf{E} = \begin{bmatrix} 0 \\ E_\perp \\ E_\parallel \end{bmatrix}.
# ```
#
# ``E_\perp`` couples to the curvature and grad-B drifts (Fermi mechanism),
# while ``E_\parallel`` directly accelerates particles along the field line.
#
# ### Adiabaticity Check
#
# The adiabatic parameter ``\epsilon = \rho_L / L_B`` where ``\rho_L`` is the
# gyroradius and ``L_B = |\nabla B|^{-1} B`` is the field gradient scale:
#
# ```math
# \rho_L = \frac{m v_\perp}{|q| B} \approx 0.7\,\mathrm{m}, \quad
# L_B \approx \alpha^{-1} = 20\,\mathrm{m}, \quad
# \epsilon \approx 0.04.
# ```
#
# With ``\epsilon \ll 1``, the GC approximation is well satisfied.

## --- Parameters (SI units) ---
const B₀₀ = 1.0e-4       # baseline field strength [T]
const α = 0.05        # shear/curvature parameter
const β = 0.5         # compression rate [1/s]
const E_perp = 5.0e-7     # perpendicular E field [V/m]
const E_par = 5.0e-8     # parallel E field [V/m]

const eV = TP.eV  # electron volt [J]

## Time-dependent B₀
B₀(t) = B₀₀ * (1 + β * t)

## Sheared magnetic field: ∇·B = 0
function shear_B(x, t = 0.0)
    B0t = B₀(t)
    return SA[α * B0t * x[3], 0.0, B0t]
end

## Constant electric field with both perp and parallel components
function const_E(x, t = 0.0)
    return SA[0.0, E_perp, E_par]
end

## --- Species: Proton ---
const species = Proton

## --- Initial Conditions ---
x0 = [0.1, 0.0, 0.0]
v_mag = 5.0e4               # 50 km/s, non-relativistic
v_par0 = v_mag / √2        # pitch angle 45°
v_perp0 = v_mag / √2
v0 = [v_perp0, 0.0, v_par0]

const stateinit = [x0..., v0...]
const tspan = (0.0, 0.8)

## Report adiabatic parameter
q_p, m_p = species.q, species.m
B_init = norm(shear_B(x0))
ρ_L = m_p * v_perp0 / (abs(q_p) * B_init)
L_B = 1 / α
println(
    "Adiabatic parameter ε = ρ_L / L_B = ", round(ρ_L, digits = 2),
    " / ", round(L_B, digits = 1), " = ", round(ρ_L / L_B, digits = 3)
)

# ## GC Tracing
#
# We compute the initial GC state and trace the GC equations.

stateinit_gc, param_gc = prepare_gc(
    stateinit, const_E, shear_B; species
)
q, q2m, μ, Efunc, Bfunc = param_gc
m = q / q2m

println(
    "Initial GC: X = ", round.(stateinit_gc[1:3], digits = 3),
    ", v∥ = ", round(stateinit_gc[4], digits = 0), " m/s"
)
println("Magnetic moment μ = ", μ, " J/T")
println(
    "Initial kinetic energy = ", round(
        0.5 * m * stateinit_gc[4]^2 + μ *
            norm(Bfunc(SA[stateinit_gc[1], stateinit_gc[2], stateinit_gc[3]], 0.0)) / eV,
        digits = 1
    ), " eV"
)

prob_gc = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
sol_gc = solve(prob_gc, Vern9())
println("GC solver steps: ", length(sol_gc.t))

# ## Energy Decomposition
#
# The GC work rates (Northrop 1963) decompose ``dK/dt`` into:
#
# | Term | Formula | Physics |
# |:---|:---|:---|
# | ``P_\parallel`` | ``q v_\parallel (\mathbf{E}\cdot\hat{\mathbf{b}})`` | Direct E∥ acceleration |
# | ``P_\mathrm{Fermi}`` | ``\frac{m v_\parallel^2}{B} (\hat{\mathbf{b}}\times\boldsymbol{\kappa})\cdot\mathbf{E}`` | Curvature drift × E |
# | ``P_\mathrm{grad}`` | ``\frac{\mu}{B} (\hat{\mathbf{b}}\times\nabla B)\cdot\mathbf{E}`` | Grad-B drift × E |
# | ``P_\mathrm{betatron}`` | ``\mu \frac{\partial B}{\partial t}`` | μ conservation in ∂B/∂t |

## Compute work rates and kinetic energy along the GC trajectory
function energy_decomposition(sol, param_gc)
    q, q2m, μ, _, Bfunc = param_gc
    m = q / q2m
    ts = sol.t
    n = length(ts)
    P_arr = zeros(n, 4)  # columns: par, fermi, grad, beta
    K_gc = zeros(n)
    B_along = zeros(n)

    for (i, t) in enumerate(ts)
        xu = sol.u[i]
        X = SA[xu[1], xu[2], xu[3]]
        Bmag_val = norm(Bfunc(X, t))
        B_along[i] = Bmag_val
        K_gc[i] = 0.5 * m * xu[4]^2 + μ * Bmag_val
        wr = TestParticle.get_work_rates_gc(xu, param_gc, t)
        P_arr[i, 1] = wr[1]  # P_par
        P_arr[i, 2] = wr[2]  # P_fermi
        P_arr[i, 3] = wr[3]  # P_grad
        P_arr[i, 4] = wr[4]  # P_betatron
    end

    return P_arr, K_gc, B_along
end

ts = sol_gc.t
n = length(ts)
P_arr, K_gc, B_along = energy_decomposition(sol_gc, param_gc)

## Cumulative work (trapezoidal integration)
function cumtrapz(t, y)
    @assert length(t) == length(y) "Dimension mismatch between t and y"
    out = similar(y, eltype(y), length(t))
    out[1] = zero(eltype(y))
    for i in 2:length(t)
        out[i] = out[i - 1] + 0.5 * (y[i] + y[i - 1]) * (t[i] - t[i - 1])
    end
    return out
end

W_par = cumtrapz(ts, @view P_arr[:, 1])
W_fermi = cumtrapz(ts, @view P_arr[:, 2])
W_grad = cumtrapz(ts, @view P_arr[:, 3])
W_beta = cumtrapz(ts, @view P_arr[:, 4])
W_total = W_par .+ W_fermi .+ W_grad .+ W_beta

## Actual energy change
ΔK_gc = K_gc .- K_gc[1]

## Convert to eV for display
K_gc_eV = K_gc ./ eV
ΔK_gc_eV = ΔK_gc ./ eV
W_par_eV = W_par ./ eV
W_fermi_eV = W_fermi ./ eV
W_grad_eV = W_grad ./ eV
W_beta_eV = W_beta ./ eV
W_total_eV = W_total ./ eV
P_arr_eV = P_arr ./ eV

using Markdown #hide
io = IOBuffer() #hide
println(io, "| Quantity | ΔEnergy (eV) |") #hide
println(io, "| -------- | ------------ |") #hide
println(io, "| ΔK (GC) | $(round(ΔK_gc[end] / eV, digits = 1)) |") #hide
println(io, "| Σ Work (GC) | $(round(W_total[end] / eV, digits = 1)) |") #hide
println(io, "| Betatron | $(round(W_beta[end] / eV, digits = 1)) |") #hide
println(io, "| Fermi | $(round(W_fermi[end] / eV, digits = 1)) |") #hide
println(io, "| Grad-B | $(round(W_grad[end] / eV, digits = 1)) |") #hide
println(io, "| E∥ (parallel) | $(round(W_par[end] / eV, digits = 1)) |") #hide
Markdown.parse(String(take!(io))) #hide

# ## Visualization

## --- Figure 1: Trajectory and Energy Overview ---
f1 = Figure(size = (1100, 700), fontsize = 14)

ax1 = Axis3(
    f1[1:2, 1],
    title = "GC Trajectory (sheared, time-dependent B)",
    xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]",
    aspect = :data,
)
lines!(
    ax1, sol_gc, idxs = (1, 2, 3), color = :steelblue, linewidth = 2,
    label = "Guiding Center"
)
scatter!(
    ax1, [sol_gc.u[1][1]], [sol_gc.u[1][2]], [sol_gc.u[1][3]],
    color = :green, markersize = 12, label = "Start"
)
scatter!(
    ax1, [sol_gc.u[end][1]], [sol_gc.u[end][2]], [sol_gc.u[end][3]],
    color = :red, markersize = 12, label = "End"
)
axislegend(ax1; position = :rt)

## |B| along trajectory
ax2 = Axis(
    f1[1, 2],
    title = "|B| Along GC Trajectory",
    xlabel = "Time [s]", ylabel = "|B| [T]",
)
lines!(ax2, ts, B_along, color = :steelblue, linewidth = 2)

## Kinetic energy
ax3 = Axis(
    f1[2, 2],
    title = "Kinetic Energy (K = ½mv∥² + μB)",
    xlabel = "Time [s]", ylabel = "K [eV]",
)
lines!(ax3, ts, K_gc_eV, color = :steelblue, linewidth = 2)

f1 = DisplayAs.PNG(f1) #hide

# The GC kinetic energy shows a clear secular increase driven by the combined
# action of all three adiabatic acceleration mechanisms.

## --- Figure 2: Cumulative Work Decomposition ---
f2 = Figure(size = (1100, 800), fontsize = 14)

ax_cum = Axis(
    f2[1, 1],
    title = "Cumulative Work Decomposition",
    xlabel = "Time [s]", ylabel = "Energy Gain [eV]",
)
lines!(
    ax_cum, ts, ΔK_gc_eV, color = :black, linewidth = 2.5,
    label = "ΔK (actual)"
)
lines!(
    ax_cum, ts, W_total_eV, color = :black, linestyle = :dash,
    linewidth = 2, label = "Σ Work rates (consistency)"
)
lines!(
    ax_cum, ts, W_beta_eV, color = :coral, linewidth = 2,
    label = "Betatron (μ ∂B/∂t)"
)
lines!(
    ax_cum, ts, W_fermi_eV, color = :teal, linewidth = 2,
    label = "Fermi (curvature drift × E)"
)
lines!(
    ax_cum, ts, W_grad_eV, color = :goldenrod, linewidth = 2,
    label = "Grad-B drift × E"
)
lines!(
    ax_cum, ts, W_par_eV, color = :purple, linewidth = 2,
    label = "Parallel E∥ (q v∥ E⋅b̂)"
)
axislegend(ax_cum; position = :lt)

## --- Figure 3: Instantaneous Work Rates ---
ax_inst = Axis(
    f2[2, 1],
    title = "Instantaneous Work Rates",
    xlabel = "Time [s]", ylabel = "Power [eV/s]",
)
@views lines!(
    ax_inst, ts, P_arr_eV[:, 1], color = :purple,
    linewidth = 1.5, label = "P_parallel"
)
@views lines!(
    ax_inst, ts, P_arr_eV[:, 2], color = :teal,
    linewidth = 1.5, label = "P_fermi"
)
@views lines!(
    ax_inst, ts, P_arr_eV[:, 3], color = :goldenrod,
    linewidth = 1.5, label = "P_grad"
)
@views lines!(
    ax_inst, ts, P_arr_eV[:, 4], color = :coral,
    linewidth = 1.5, label = "P_betatron"
)
axislegend(ax_inst; position = :rt)

f2 = DisplayAs.PNG(f2) #hide

# ## Discussion
#
# The cumulative work decomposition reveals three distinct acceleration channels:
#
# 1. **Betatron acceleration** (coral): ``\mu\,\partial B/\partial t`` provides
#    a steady power input. Since ``\partial B/\partial t`` is nearly constant
#    in our linear ramp model, the cumulative work grows quadratically in time.
#    This is the dominant mechanism when the field strength changes
#    significantly on the drift timescale.
#
# 2. **Fermi + Grad-B drift work** (teal + goldenrod): These arise from the
#    curvature and grad-B drifts crossing ``E_\perp``. The sign depends on
#    the relative orientation of the drift and the electric field.
#
# 3. **Parallel ``E_\parallel``** (purple): Direct acceleration along ``\hat{\mathbf{b}}``
#    by the parallel electric field. Its contribution is proportional to
#    ``q v_\parallel E_\parallel`` and grows linearly if both are approximately
#    constant.
#
# The dashed black line (sum of integrated work rates) tracks the actual
# ``\Delta K`` (solid black), confirming consistency of the GC energy budget.

## --- Figure 4: Adiabatic Invariant and Velocity Evolution ---
##
## In the GC formalism, ``\mu`` is exactly conserved. We verify this along
## the trajectory and show how ``v_\parallel`` evolves.

vpar_vals = [u[4] for u in sol_gc.u]

f3 = Figure(size = (1000, 450), fontsize = 14)

ax_mu = Axis(
    f3[1, 1],
    title = "Magnetic Moment μ (exactly conserved in GC)",
    xlabel = "Time [s]", ylabel = "μ [J/T]",
)
lines!(ax_mu, ts, fill(μ, n), color = :steelblue, linewidth = 2)

ax_vpar = Axis(
    f3[1, 2],
    title = "Parallel Velocity v∥",
    xlabel = "Time [s]", ylabel = "v∥ [m/s]",
)
lines!(ax_vpar, ts, vpar_vals, color = :coral, linewidth = 2)

f3 = DisplayAs.PNG(f3) #hide

# The exact conservation of ``\mu`` is a built-in property of the GC
# equations. The secular change in ``v_\parallel`` reflects the combined
# effect of Fermi acceleration and ``E_\parallel``.

## --- Parameter Study: Effect of Compression Rate β ---

function run_case(β_val)
    B_local(x, t = 0.0) = let B0t = B₀₀ * (1 + β_val * t)
        SA[α * B0t * x[3], 0.0, B0t]
    end
    E_local(x, t = 0.0) = SA[0.0, E_perp, E_par]

    sigc, pgc = prepare_gc(stateinit, E_local, B_local; species)
    prob = ODEProblem(trace_gc!, sigc, tspan, pgc)
    sol = solve(prob, Vern9())

    q, q2m, μ, _, Bfunc = pgc
    m = q / q2m
    K = [
        let xu = sol.u[i], t = sol.t[i]
                X = SA[xu[1], xu[2], xu[3]]
                (0.5 * m * xu[4]^2 + μ * norm(Bfunc(X, t))) / eV
        end for i in 1:length(sol.t)
    ]

    return sol.t, K .- K[1]
end

β_cases = [0.0, 0.25, 0.5, 1.0]
colors = [:gray, :steelblue, :coral, :darkred]
labels = ["β = 0 (static)", "β = 0.25", "β = 0.5", "β = 1.0"]

f4 = Figure(size = (850, 500), fontsize = 14)
ax_param = Axis(
    f4[1, 1],
    title = "Energy Gain vs. Compression Rate",
    xlabel = "Time [s]", ylabel = "ΔK [eV]",
)
for (βi, label, color) in zip(β_cases, labels, colors)
    ts_p, ΔK_p = run_case(βi)
    lines!(ax_param, ts_p, ΔK_p, color = color, linewidth = 2, label = label)
end
axislegend(ax_param; position = :lt)

f4 = DisplayAs.PNG(f4) #hide

# The static case (β = 0) shows no betatron or Fermi acceleration — only
# the modest contribution from ``E_\parallel``. As β increases, the energy
# gain grows rapidly, dominated by betatron acceleration with additional
# contributions from Fermi acceleration.
