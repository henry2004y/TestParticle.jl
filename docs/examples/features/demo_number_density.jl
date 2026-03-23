# # Detector Flux Calculation
#
# This example demonstrates how to calculate the particle flux across multiple detectors
# (planes) and compares the results with the analytical solution for a freely expanding
# Maxwellian cloud.

import DisplayAs #hide
using TestParticle
import TestParticle as TP
using VelocityDistributionFunctions
using OrdinaryDiffEq
using StaticArrays
using Statistics
using LinearAlgebra
using Random
using CairoMakie
using Meshes
CairoMakie.activate!(type = "png") #hide

Random.seed!(1234);

# ## Simulation Setup

# Parameters
N = 100_000
m = 1.0
q = 1.0
T = 1.0
# Thermal velocity: ``v_{th} = \sqrt{2 k_B T/m}``
vth = sqrt(2 * T / m)
t_end = 2.0

# Initialize particles
# Point source at origin: ``\mathbf{r}_0 = \mathbf{0}``
x0 = [SVector(0.0, 0.0, 0.0) for _ in 1:N];

# Maxwellian velocity distribution
# ``u_{bulk} = 0``, ``n=1`` effectively for distribution shape
vdf = TP.Maxwellian([0.0, 0.0, 0.0], T, 1.0; m = m)
# Use the vectorized rand to get a Vector of SVectors
v0 = rand(vdf, N);

# Create `TraceProblem` template
# Construct an `EnsembleProblem` to simulate multiple trajectories.
prob_func(prob, i, repeat) = remake(prob, u0 = vcat(x0[i], v0[i]))

# Define a single problem template
u0_dummy = vcat(x0[1], v0[1])
## Zero fields
param = prepare(TP.ZeroField(), TP.ZeroField(); q, m)
tspan = (0.0, t_end)
prob = ODEProblem(trace, u0_dummy, tspan, param);

ensemble_prob = EnsembleProblem(prob; prob_func)

# ## Solve
sols = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories = N);

# ## Flux Calculation with Multiple Detectors

# We define a set of planes along the x-axis to detect particle crossings.
L = 6.0
num_planes = 40
xs = range(0.1, L, length = num_planes)
detectors = [Plane(Point(x, 0.0, 0.0), Vec(1.0, 0.0, 0.0)) for x in xs]

# Calculate flux across all planes
# We use a weight of 1.0 per test particle.
weights = 1.0
velocity_fluxes, n_fluxes = get_particle_fluxes(sols, detectors, weights);

# ## Analytical Solution
#
# For a point source at origin at ``t=0``, a particle with velocity ``v_x``
# will cross the plane at ``x_c > 0`` at time ``t = x_c / v_x``.
# The particle crosses the plane before ``t_{end}`` if ``x_c / v_x < t_{end}``,
# which means ``v_x > x_c / t_{end}``.
#
# The fraction of particles with ``v_x > v_{xc}`` in a Maxwellian distribution is:
# ``P(v_x > v_{xc}) = \frac{1}{\sqrt{\pi} v_{th}} \int_{v_{xc}}^{\infty} \exp(-v_x^2/v_{th}^2) dv_x = \frac{1}{2} \text{erfc}(v_{xc}/v_{th})``
#
# Thus, the expected number of particles crossing the plane at ``x_c`` is:
# ``N_{cross} = \frac{N}{2} \text{erfc}\left(\frac{x_c}{v_{th} t_{end}}\right)``

using SpecialFunctions: erfc
n_analytic = [N / 2 * erfc(x / (vth * t_end)) for x in xs]

# ## Visualization

f = Figure(fontsize = 18)
ax = Axis(
    f[1, 1],
    title = "Cumulative Particle Crossings (t ≤ $(t_end))",
    xlabel = "Plane Position (x)",
    ylabel = "Cumulative Counts"
)

scatter!(ax, xs, n_fluxes, label = "Simulation", color = :blue, markersize = 10)
lines!(ax, xs, n_analytic, label = "Analytical", color = :red, linewidth = 2)
axislegend(ax)

f = DisplayAs.PNG(f) #hide
