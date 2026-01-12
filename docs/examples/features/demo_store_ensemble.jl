# # Storage of Ensemble Solutions
#
# When running large ensemble simulations, storing the full high-precision solution can
# consume a significant amount of memory and disk space. This example demonstrates two
# strategies to reduce storage requirements:
#
# 1. Using `strip_solution` to remove dense interpolation data before saving.
# 2. Converting trajectory data to `Float32` before saving.
#
# We will use `JLD2.jl` for saving the data.

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using SciMLBase
using StaticArrays
using JLD2
using Random

Random.seed!(1234);

# ## Setup: Ensemble Simulation
#
# First, we set up a simple ensemble problem: protons in a uniform magnetic field.

B_uniform(x) = SA[0, 0, 1.0e-8]
E_zero = TestParticle.ZeroField()

## Initialize a proton with some perpendicular velocity
x0 = [0.0, 0.0, 0.0]
u0 = [1.0, 0.0, 0.0]
stateinit = [x0..., u0...]

param = prepare(E_zero, B_uniform, species = Proton)
tspan = (0.0, 10.0)

## Define the problem
prob = ODEProblem(trace!, stateinit, tspan, param)

## Define a prob_func to vary initial velocity slightly
prob_func(prob, i, repeat) = remake(prob, u0 = prob.u0 .* (1 + 0.1 * rand()))

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)

## Solve for a small ensemble
## We use a high-order solver (Vern9) which produces dense interpolation data by default
trajectories = 2
sols = solve(ensemble_prob, Vern9(), EnsembleThreads(); trajectories, saveat = 0.2);

# ## Case 1: Converting to Float32
#
# For visualization or statistical analysis, `Float64` precision is often unnecessary.
# We can convert the position and time data to `Float32` before saving.
#
# Note: Since `sols` is an `EnsembleSolution` containing a vector of `ODESolution`s,
# we iterate through them.

filename_f32 = tempname() * ".jld2"

## Extract and convert data
## We structure the data as a vector of named tuples or a dictionary for saving
output_data = map(sols) do s
    (t = Float32.(s.t), u = [Float32.(state) for state in s.u])
end

## Save to JLD2
jldsave(filename_f32; data = output_data)

## Verify loading
loaded_f32 = load(filename_f32)["data"]
println("Loaded Float32 data type: ", eltype(loaded_f32[1].t))
println("Float32 conversion successfully saved and loaded.")

# ## Cleanup
#
# Remove the temporary files.

rm(filename_f32, force = true)
