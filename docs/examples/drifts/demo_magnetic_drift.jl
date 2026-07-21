# # Magnetic Drift and Energy Partition
#
# This example demonstrates the combined magnetic drift (curvature + grad-B) of a
# single charged particle and how it partitions between the parallel and
# perpendicular kinetic energies. It is the direct numerical test of the question
# raised in [issue #61](https://github.com/henry2004y/TestParticle.jl/issues/61):
# in a field where curvature and grad-B are coupled, are the two drifts a single
# "magnetic drift" distinguished only by the energy partition?
#
# More theoretical details can be found in
# [Magnetic Drift](https://henry2004y.github.io/KeyNotes/contents/single.html#b-b-grad-b-drift)
# and *Fundamentals of Plasma Physics* by Paul Bellan, §3.6.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
using LinearAlgebra: norm, ×, ⋅, normalize
using ForwardDiff: jacobian
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Field and drift deduction
#
# We use a 2D azimuthal field that circles the z-axis with radius `r`:
# `B = B₀ (x₂, -x₁, 0) / r²`. It has both a curved field line (curvature radius
# `R_c = r`) and a perpendicular magnitude gradient (`|B| = B₀/r`), and the two are
# coupled through `B·R_c = B₀ = const`. With `E = 0` the guiding-center velocity
# (mirroring [`trace_gc_drifts!`](@ref)) is
# ```math
# \mathbf{v}_\mathrm{gc} =
#   \frac{m}{2qB^2}(\mathbf{b}\times\nabla B)\,v_\perp^2 +
#   \frac{m}{qB}(\mathbf{b}\times\boldsymbol{\kappa})\,v_\parallel^2 ,
# \qquad \boldsymbol{\kappa} = (\mathbf{b}\cdot\nabla)\mathbf{b}.
# ```
# For this field both `(b×∇B)` and `(b×κ)` point along `-ẑ` and their denominators
# collapse to the same constant `B₀`, giving the clean combined form
# ```math
# v_z = -\frac{m}{qB_0}\left(v_\parallel^2 + \frac{v_\perp^2}{2}\right).
# ```
# So curvature drift carries the `v_∥²` part and grad-B drift carries the
# `v_⊥²/2` part of one and the same magnetic drift.

const Bmag0 = 1.0e-7

curve_B(x) = SA[x[2] / norm(x[1:2])^2, -x[1] / norm(x[1:2])^2, 0.0] * Bmag0
zero_E(x) = SA[0.0, 0.0, 0.0]

# Analytic guiding-center drift evaluated directly from the field function, using
# exactly the same quantities (`∇B` from the B-Jacobian, `κ` from `JB·b`) as the
# package's [`get_B_parameters`](@ref) / `trace_gc_drifts!`. This is the "code
# tracing" reference against which we compare the Boris-traced orbit.
function analytic_drift(x, v, q, m)
    B = curve_B(x)
    Bmag = norm(B)
    b = B / Bmag
    JB = jacobian(curve_B, x)
    ∇B = JB' * b
    κ = (JB * b + b * (-(∇B ⋅ b))) / Bmag
    Ω = (q / m) * Bmag
    vpar = v ⋅ b
    vperp = v - vpar * b
    wsq = vperp ⋅ vperp
    return wsq * (b × ∇B) / (2Ω * Bmag) + vpar^2 * (b × κ) / Ω
end

# ## Pitch-angle scan
#
# We launch protons from the same point `x0` with a fixed speed `v0` but varying
# pitch angle `α` between `v` and `b`, i.e. `v_∥ = v0 cos α`, `v_⊥ = v0 sin α`.
# The predicted drift is then compared to the drift of the guiding center obtained
# from the full Boris orbit.

x0 = SA[1.0, 0.0, 0.0]
b0 = normalize(curve_B(x0))          # b = (0, -1, 0) at x0
eperp = SA[0.0, 0.0, 1.0]             # perpendicular to b at x0
const v0 = 1.0
tspan = (0.0, 30.0)
n = 9
αs = range(0, π / 2; length = n)

q, m = Proton.q, Proton.m
vpar_list = Float64[]
vperp_list = Float64[]
theory_z = Float64[]
meas_z = Float64[]
sols = []   # keep solutions for the trajectory plot

param = prepare(zero_E, curve_B, species = Proton)
gc = get_gc_func(param)
ts = range(tspan..., length = 400)
A = hcat(ones(length(ts)), ts)

for α in αs
    vpar = v0 * cos(α)
    vperp = v0 * sin(α)
    v = vpar * b0 + vperp * eperp
    push!(vpar_list, vpar)
    push!(vperp_list, vperp)

    # Theoretical prediction from the field itself.
    push!(theory_z, analytic_drift(x0, v, q, m)[3])

    # Boris trace.
    stateinit = [x0..., v...]
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Vern9(); abstol = 1.0e-10, reltol = 1.0e-10)
    push!(sols, sol)

    # Drift of the guiding center: slope of its z position vs time.
    z_gc = [gc(sol(t))[3] for t in ts]
    push!(meas_z, (A \ z_gc)[2])
end

# ## Results
#
# Scatter of the Boris-traced guiding-center drift against the analytic drift
# computed directly from the field function (should lie on the `y = x` line), and
# the linear dependence on the energy partition `(v_∥² + v_⊥²/2)`.

fig = Figure(size = (1500, 500), fontsize = 18)

ax3 = Axis3(
    fig[1, 1]; title = "Sample trajectories", xlabel = "x", ylabel = "y",
    zlabel = "z", aspect = :data
)
for idx in (1, 5, 9)   # α = 0 (curvature only), 45° (mixed), 90° (grad-B only)
    lines!(
        ax3, sols[idx]; idxs = (1, 2, 3),
        label = "α = $(round(rad2deg(αs[idx]); digits = 0))°"
    )
end
axislegend(ax3, framevisible = true, backgroundcolor = (:white, 0.6))

ax_conf = Axis(
    fig[1, 2]; title = "Code tracing confirmation",
    xlabel = "Analytic drift vz [m/s]", ylabel = "Measured drift vz [m/s]",
    aspect = 1
)
scatter!(ax_conf, theory_z, meas_z; markersize = 10)
xr = extrema(theory_z)
lines!(
    ax_conf, [xr[1], xr[2]], [xr[1], xr[2]]; color = :red, linestyle = :dash,
    label = "y = x"
)
axislegend(ax_conf, position = :lt)

ax_part = Axis(
    fig[1, 3]; title = "Energy partition",
    xlabel = "v∥² + v⊥²/2  [m²/s²]", ylabel = "Drift vz [m/s]"
)
X = vpar_list .^ 2 .+ 0.5 .* (vperp_list .^ 2)
scatter!(ax_part, X, meas_z; markersize = 10, color = :black, label = "measured")
slope = (X' * X) \ (X' * meas_z)   # least-squares slope
xs = range(extrema(X)...; length = 50)
lines!(ax_part, xs, slope[1] .* xs; color = :blue, label = "linear fit")
axislegend(ax_part, position = :lt)

# The fitted slope equals m/(q B₀) to within the (small) finite-Larmor-radius
# correction; the grad-B and curvature parts therefore carry the v⊥²/2 and v∥²
# weights, i.e. they are one magnetic drift, not two independent forces.
fig = DisplayAs.PNG(fig) #hide
