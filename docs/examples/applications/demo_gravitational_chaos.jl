# # Gravitational Chaos and Fractal Basins
#
# Three equal masses can follow a circular equilateral-orbit solution of the
# Newtonian three-body problem. We use their time-dependent gravity to trace
# massless test particles. Nearby initial conditions can separate rapidly, and
# the final nearest-primary map reveals intricate basin boundaries.

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png") #hide

const G = 1.0
const primary_mass = 1.0
const orbit_radius = 1.0
const softening = 0.05
const angular_frequency = sqrt(
    3G * primary_mass / (3orbit_radius^2 + softening^2)^(3 / 2)
)
const orbital_period = 2π / angular_frequency
const final_time = 3orbital_period

function primary_position(i, t)
    phase = angular_frequency * t + 2π * (i - 1) / 3
    return SA[
        orbit_radius * cos(phase),
        orbit_radius * sin(phase),
        0.0
    ]
end

function gravity(x, t)
    r = SA[x[1], x[2], x[3]]
    acceleration = zero(r)
    for i in 1:3
        displacement = r - primary_position(i, t)
        distance_squared = displacement ⋅ displacement + softening^2
        acceleration -= G * primary_mass * displacement / distance_squared^(3 / 2)
    end
    return acceleration
end

param = prepare(ZeroField(), ZeroField(), gravity; q = 0.0, m = 1.0)

# ## Fractal-like outcome map
#
# Each grid point is an initial position for a particle at rest. Its color
# indicates which primary it is closest to after three orbital periods. The
# finite-time map approximates the basins; increasing the grid resolution and
# integration time reveals finer structure.

resolution = 32
xgrid = range(-2.0, 2.0; length = resolution)
ygrid = range(-2.0, 2.0; length = resolution)
initial_states = [
    SA[x, y, 0.0, 0.0, 0.0, 0.0]
    for y in ygrid for x in xgrid
]
tspan = (0.0, final_time)
base_prob = ODEProblem(trace, first(initial_states), tspan, param)

prob_func = (prob, ctx) -> remake(prob; u0 = initial_states[ctx.sim_id])
output_func = (sol, i) -> (sol.u[end], false)
ensemble_prob = EnsembleProblem(
    base_prob; prob_func, output_func, safetycopy = false
)
endpoints = solve(
    ensemble_prob, Tsit5(), EnsembleSerial();
    trajectories = length(initial_states),
    save_everystep = false, save_start = false, save_end = true
)

outcomes = Matrix{Int}(undef, resolution, resolution)
for (index, endpoint) in enumerate(endpoints.u)
    ix = mod1(index, resolution)
    iy = div(index - 1, resolution) + 1
    final_position = SA[endpoint[1], endpoint[2], endpoint[3]]
    distances = [
        norm(final_position - primary_position(i, final_time))
        for i in 1:3
    ]
    outcomes[ix, iy] = argmin(distances)
end

# ## Trajectory sensitivity
#
# The next pair starts with positions separated by only `1e-5` in the x
# direction. Their different paths illustrate sensitivity to initial conditions.

nearby_states = (
    SA[1.45, 0.0, 0.0, 0.0, 0.25, 0.0],
    SA[1.45 + 1.0e-5, 0.0, 0.0, 0.0, 0.25, 0.0]
)
sample_times = range(0.0, final_time; length = 1500)
sample_solutions = [
    solve(
        remake(base_prob; u0), Tsit5();
        saveat = sample_times, abstol = 1.0e-9, reltol = 1.0e-9
    )
    for u0 in nearby_states
]

# ## Visualization

fig = Figure(size = (1200, 560))
ax_map = Axis(
    fig[1, 1];
    title = "Finite-time nearest-primary basins",
    xlabel = "Initial x", ylabel = "Initial y", aspect = DataAspect()
)
basins = heatmap!(
    ax_map, xgrid, ygrid, outcomes;
    colormap = [:dodgerblue, :darkorange, :seagreen],
    colorrange = (1, 3), interpolate = false
)
Colorbar(
    fig[1, 2], basins;
    ticks = (1:3, ["Primary 1", "Primary 2", "Primary 3"]),
    label = "Nearest primary"
)

ax_orbit = Axis(
    fig[1, 3];
    title = "Nearby test-particle trajectories",
    xlabel = "x", ylabel = "y", aspect = DataAspect()
)
primary_times = range(0.0, orbital_period; length = 300)
for i in 1:3
    orbit = [primary_position(i, t) for t in primary_times]
    lines!(
        ax_orbit, getindex.(orbit, 1), getindex.(orbit, 2);
        color = Makie.wong_colors()[i], linestyle = :dash, label = "Primary $i"
    )
end
for (i, sol) in enumerate(sample_solutions)
    lines!(
        ax_orbit, sol[1, :], sol[2, :];
        linewidth = 2, label = "Test particle $i"
    )
end
axislegend(ax_orbit; position = :rt, labelsize = 11)

fig = DisplayAs.PNG(fig) #hide

# This is a softened, dimensionless model. The outcome map is finite-time and
# resolution-dependent; its detailed boundaries should not be interpreted as
# an exact infinite-time fractal.
