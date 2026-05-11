# # Polarization Drift
#
# This example demonstrates a single proton motion under time-varying E field.
# More theoretical details can be found in [Time-Varying E Drift](https://henry2004y.github.io/KeyNotes/contents/single.html#sec-time-varying_E), and Fundamentals of Plasma Physics by Paul Bellan.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Plotting Function
#
# We define a function to visualize the results.

function plot_drift_case(sol, title)
    fig = Figure(size = (1200, 600), fontsize = 20)

    ## 1. Left Column: Time Series
    gl_left = fig[1, 1] = GridLayout()
    ax_pos = Axis(gl_left[1, 1], ylabel = "Position [m]", title = title)
    lines!(ax_pos, sol, idxs = (0, 1), label = "x")
    lines!(ax_pos, sol, idxs = (0, 2), label = "y")
    lines!(ax_pos, sol, idxs = (0, 3), label = "z")
    axislegend(ax_pos, position = :lt, framevisible = true, backgroundcolor = (:white, 0.5))
    hidexdecorations!(ax_pos, grid = false)

    ax_vel = Axis(gl_left[2, 1], xlabel = "Time [s]", ylabel = "Velocity [m/s]")
    lines!(ax_vel, sol, idxs = (0, 4), label = "vx")
    lines!(ax_vel, sol, idxs = (0, 5), label = "vy")
    lines!(ax_vel, sol, idxs = (0, 6), label = "vz")
    axislegend(ax_vel, position = :lt, framevisible = true, backgroundcolor = (:white, 0.5))

    linkxaxes!(ax_pos, ax_vel)

    ## 2. Right Column: 3D Trajectory
    ax_3d = Axis3(
        fig[1, 2];
        title = "3D Trajectory", xlabel = "x", ylabel = "y", zlabel = "z", aspect = :data
    )

    gc = sol.prob.p[1] |> p -> get_gc_func(p; use_vE = true)
    gc_plot(x, y, z, vx, vy, vz, t) = (gc(SA[x, y, z, vx, vy, vz, t])...,)

    lines!(
        ax_3d, sol;
        idxs = (1, 2, 3), color = Makie.wong_colors()[1], alpha = 0.5, label = "Particle"
    )
    lines!(
        ax_3d, sol;
        idxs = (gc_plot, 1, 2, 3, 4, 5, 6, 0),
        color = Makie.wong_colors()[2], label = "GC from Orbit"
    )
    axislegend(ax_3d, framevisible = true, backgroundcolor = (:white, 0.5))

    return fig
end

# ## Time-Varying E Field
#
# We use a time-varying electric field and a uniform magnetic field and trace a proton.

uniform_B(x) = SA[0, 0, 1.0e-8]

function time_varying_E(x, t)
    return SA[0, 1.0e-9 * 0.1 * t, 0]
end

## Initial condition
stateinit = let x0 = [1.0, 0, 0], v0 = [0.0, 1.0, 0.1]
    [x0..., v0...]
end
## Time span
tspan = (0, 100)
param = prepare(time_varying_E, uniform_B, species = Proton)
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern9())

fig = plot_drift_case(sol, "Polarization Drift Case")
fig = DisplayAs.PNG(fig) #hide
