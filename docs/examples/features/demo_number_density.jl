# # Number Density Calculation
#
# This example demonstrates how to calculate the number density of a group of particles
# and compares it with the analytical solution for a freely expanding Maxwellian cloud.

using TestParticle
using Meshes
using OrdinaryDiffEq
using StaticArrays
using Statistics
using LinearAlgebra
using Random
using CairoMakie
CairoMakie.activate!(type = "png") #hide
import DisplayAs #hide

Random.seed!(1234)

# ## Simulation Setup

# Parameters
N = 100_000
m = 1.0
q = 1.0
T = 1.0
# Thermal velocity in TestParticle (and VelocityDistributionFunctions) is vth = sqrt(2*T/m)
vth = sqrt(2 * T / m)
t_end = 2.0

# Initialize particles
# Point source at origin
x0 = [SVector(0.0, 0.0, 0.0) for _ in 1:N]

# Maxwellian velocity distribution
# u0 = 0, p = n*T = 1*1 = 1 (assuming n=1 effectively for distribution shape), n=1
vdf = TestParticle.Maxwellian([0.0, 0.0, 0.0], T, 1.0; m=m)
# Use the vectorized rand to get a Vector of SVectors
v0 = rand(vdf, N)

# Create TraceProblem template
# We need to create a prob for each particle or use an ensemble.
# TestParticle.trace usually takes one particle.
# We will construct an EnsembleProblem.

function prob_func(prob, i, repeat)
    remake(prob, u0 = vcat(x0[i], v0[i]))
end

# Define a single problem template
u0_dummy = vcat(x0[1], v0[1])
# Zero fields
param = prepare(TestParticle.ZeroField(), TestParticle.ZeroField(), species=User, q=q, m=m)
prob = TraceProblem(u0_dummy, (0.0, t_end), param)

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

# Solve
# We use Vern9 for high accuracy, though free expansion is exact.
sols = solve(ensemble_prob, Vern9(), EnsembleThreads(), trajectories=N);

# ## Density Calculation

# Define a grid for density calculation
# The cloud expands to r ~ vth * t_end.
# vth = sqrt(2) ~ 1.414. r ~ 1.414 * 2 ~ 2.8.
# Let's cover [-6, 6] to see the tails.
L = 6.0
dims = (50, 50, 50)
grid = CartesianGrid((-L, -L, -L), (L, L, L), dims=dims)

# Calculate density at t_end
# The function expects a vector of solutions.
# get_number_density returns count / volume
density = TestParticle.get_number_density(sols, grid, t_end)

# Extract a 1D slice along x-axis (y=0, z=0 approx)
mid_y = dims[2] ÷ 2
mid_z = dims[3] ÷ 2
density_x = density[:, mid_y, mid_z]

# Grid coordinates for plotting
# makegrid returns ranges for each dimension
grid_x, grid_y, grid_z = TestParticle.makegrid(grid)
xs_plot = collect(grid_x)

# ## Analytical Solution

# For a point source expanding with Maxwellian velocities:
# The spatial distribution at time t is a Gaussian with standard deviation sigma_x = vth * t / sqrt(2).
# Actually, let's derive carefully from the VDF.
# f(v) = (1 / (pi * vth^2)^(3/2)) * exp(-v^2 / vth^2)
# r = v * t => v = r / t.
# Jacobian from v-space to r-space is |dv/dr| = 1/t^3.
# n(r, t) = N * f(r/t) * (1/t^3)
#         = N * (1 / (pi * vth^2)^(3/2)) * exp(-(r/t)^2 / vth^2) / t^3
#         = N / (pi * (vth*t)^2)^(3/2) * exp(-r^2 / (vth*t)^2)
#
# Note: The Gaussian form is exp(-r^2 / (2*sigma^2)).
# Here we have exp(-r^2 / (vth*t)^2).
# So 2*sigma^2 = (vth*t)^2 => sigma = vth*t / sqrt(2).
#
# Let sigma_param = vth * t_end. The exponential is exp(-r^2 / sigma_param^2).

r2 = xs_plot.^2
# We are taking a slice at y~0, z~0. So r^2 = x^2.
# Note: The grid cells have finite volume. get_number_density returns density.
# So we compare directly with n(x, 0, 0, t).

sigma_param = vth * t_end
n_analytic = @. N / (sqrt(π) * sigma_param)^3 * exp(-r2 / sigma_param^2)

# ## Visualization

f = Figure(fontsize=18)
ax = Axis(f[1, 1],
    title = "Free Expansion Density (t = $(t_end))",
    xlabel = "x",
    ylabel = "Number Density")

scatter!(ax, xs_plot, density_x, label="Simulation", color=:blue, markersize=10)
lines!(ax, xs_plot, n_analytic, label="Analytical", color=:red, linewidth=2)
axislegend(ax)

f = DisplayAs.PNG(f) #hide
