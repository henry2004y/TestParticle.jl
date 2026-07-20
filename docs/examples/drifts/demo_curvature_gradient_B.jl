# # Curvature and Grad-B Drifts
#
# This example demonstrates a single proton motion under two different
# non-uniform magnetic fields — one with field-line curvature, one with a
# perpendicular magnitude gradient. Both are more complex than the simpler
# [ExB Drift](@ref). For each field we trace the full orbit and compare the
# guiding center (GC) obtained from the orbit with the analytic GC.
#
# - **Curvature drift** comes from a curved field line, `κ = (b·∇)b`; the GC
#   is traced with [`trace_gc_drifts!`](@ref).
# - **Grad-B drift** comes from `∇B ⊥ B`; here the analytic GC includes the
#   ExB drift too. Note that the orbit-derived and analytic GCs differ, because
#   the GC has higher-order terms beyond the first-order textbook formula.
#
# The combined effect of the two drifts for one particle is validated in
# [`demo_magnetic_drift`](@ref), and their relation at the particle and fluid
# level is discussed in [`demo_drift_relation`](@ref).
#
# More theoretical details can be found in
# [Curvature/Grad-B Drift](https://henry2004y.github.io/KeyNotes/contents/single.html#b-b-grad-b-drift),
# and *Fundamentals of Plasma Physics* by Paul Bellan.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
using LinearAlgebra: ×, ⋅, normalize, norm
using ForwardDiff: gradient
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Plotting function
#
# For a traced orbit `sol` and its analytic GC `sol_gc`, draw the time series on
# the left and the 3D trajectory (particle, GC from orbit, analytic GC) on the
# right, into the given figure row.

function plot_drift_case(fig, row, sol, sol_gc, title)
    gl_left = fig[row, 1] = GridLayout()
    ax_pos = Axis(gl_left[1, 1], ylabel = "Position [m]", title = title)
    lines!(ax_pos, sol, idxs = (0, 1), label = "x")
    lines!(ax_pos, sol, idxs = (0, 2), label = "y")
    lines!(ax_pos, sol, idxs = (0, 3), label = "z")
    axislegend(
        ax_pos, position = :lt, framevisible = true,
        backgroundcolor = (:white, 0.5)
    )
    hidexdecorations!(ax_pos, grid = false)

    ax_vel = Axis(gl_left[2, 1], xlabel = "Time [s]", ylabel = "Velocity [m/s]")
    lines!(ax_vel, sol, idxs = (0, 4), label = "vx")
    lines!(ax_vel, sol, idxs = (0, 5), label = "vy")
    lines!(ax_vel, sol, idxs = (0, 6), label = "vz")
    axislegend(
        ax_vel, position = :lt, framevisible = true,
        backgroundcolor = (:white, 0.5)
    )

    linkxaxes!(ax_pos, ax_vel)

    ax_3d = Axis3(
        fig[row, 2]; title = "3D Trajectory", xlabel = "x",
        ylabel = "y", zlabel = "z", aspect = :data
    )

    gc = get_gc_func(sol_gc.prob.p)
    gc_plot(x, y, z, vx, vy, vz) = (gc(SA[x, y, z, vx, vy, vz])...,)

    lines!(
        ax_3d, sol; idxs = (1, 2, 3), color = Makie.wong_colors()[1],
        alpha = 0.5, label = "Particle"
    )
    lines!(
        ax_3d, sol; idxs = (gc_plot, 1, 2, 3, 4, 5, 6),
        color = Makie.wong_colors()[2], label = "GC from Orbit"
    )
    lines!(
        ax_3d, sol_gc; idxs = (1, 2, 3), color = Makie.wong_colors()[3],
        linewidth = 3, label = "Analytic GC"
    )
    axislegend(ax_3d, framevisible = true, backgroundcolor = (:white, 0.5))

    return nothing
end

# ## Curvature B field
#
# A 2D azimuthal field circling the z-axis with radius `r`:
# `B = B₀ (x₂, -x₁, 0) / r²`. It has a curved field line (`R_c = r`) and we
# trace a proton with a small `v_z` so curvature drift is visible.

curved_B(x) = SA[x[2] / norm(x[1:2])^2, -x[1] / norm(x[1:2])^2, 0.0] * 1.0e-8
zero_E(x) = SA[0.0, 0.0, 0.0]

stateinit_curv = let x0 = [1.0, 0.0, 0.0], v0 = [0.0, 1.0, 0.1]
    [x0..., v0...]
end
tspan_curv = (0, 40)
param_curv = prepare(zero_E, curved_B, species = Proton)
prob_curv = ODEProblem(trace!, stateinit_curv, tspan_curv, param_curv)
sol_curv = solve(prob_curv, Vern9())
gc_curv = param_curv |> get_gc_func
gc_x0_curv = gc_curv(stateinit_curv) |> Vector
prob_gc_curv = ODEProblem(
    trace_gc_drifts!, gc_x0_curv, tspan_curv,
    (param_curv..., sol_curv)
)
sol_gc_curv = solve(prob_gc_curv, Vern7(); save_idxs = [1, 2, 3])

# ## Grad-B field
#
# A straight field `B = B₀(1 + g·x₂) ẑ` with a perpendicular gradient, plus a
# small uniform `E` so the analytic GC also carries the ExB drift. The GC is
# traced with a hand-written drift law that keeps the first-order grad-B and ExB
# terms.

grad_B(x) = SA[0, 0, 1.0e-8 + 1.0e-9 * x[2]]
uniform_E(x) = SA[1.0e-9, 0, 0]
abs_B(x) = norm(grad_B(x))

function trace_gc!(dx, x, p, t)
    q2m, _, E, B, _, sol = p
    xu = sol(t)
    gradient_B = gradient(abs_B, x)
    Bv = B(x)
    b = normalize(Bv)
    v_par = (xu[4:6] ⋅ b) .* b
    v_perp = xu[4:6] - v_par
    return dx[1:3] = norm(v_perp)^2 * (Bv × gradient_B) / (2 * q2m * norm(Bv)^3) +
        (E(x) × Bv) / norm(Bv)^2 + v_par
end

stateinit_grad = let x0 = [1.0, 0, 0], v0 = [0.0, 1.0, 0.1]
    [x0..., v0...]
end
tspan_grad = (0, 20)
param_grad = prepare(uniform_E, grad_B, species = Proton)
prob_grad = ODEProblem(trace!, stateinit_grad, tspan_grad, param_grad)
sol_grad = solve(prob_grad, Vern9())
gc_grad = param_grad |> get_gc_func
gc_x0_grad = gc_grad(stateinit_grad) |> Vector
prob_gc_grad = ODEProblem(
    trace_gc!, gc_x0_grad, tspan_grad,
    (param_grad..., sol_grad)
)
sol_gc_grad = solve(prob_gc_grad, Vern7(); save_idxs = [1, 2, 3])

# ## Results
#
# Top: curvature drift. Bottom: grad-B drift. In the grad-B case the orbit-derived
# and analytic GCs follow different trajectories, because the GC carries
# higher-order terms beyond the first-order drift formula.

fig = Figure(size = (1200, 1100), fontsize = 18)
plot_drift_case(fig, 1, sol_curv, sol_gc_curv, "Curvature Drift Case")
plot_drift_case(fig, 2, sol_grad, sol_gc_grad, "Grad-B Drift Case")
fig = DisplayAs.PNG(fig) #hide
