using TestParticle
using StaticArrays
using BenchmarkTools
using OrdinaryDiffEq

# Define a simple time-independent field
function B_field(x, t)
    return SVector(0.0, 0.0, 1e-9)
end

function E_field(x, t)
    return SVector(0.0, 0.0, 0.0)
end

# Define a time-dependent field
function B_field_td(x, t)
    return SVector(0.0, 0.0, 1e-9 * cos(t))
end

# Setup problem
x0 = [1.0, 0.0, 0.0] # Initial position
v0 = [0.0, 1.0, 0.0] # Initial velocity
u0 = [x0; v0]
tspan = (0.0, 1.0)

# Prepare fields
# Using functions directly
param = prepare(E_field, B_field; species=Proton)
prob = TraceProblem(u0, tspan, param)

println("Benchmarking time-independent field with current implementation:")
@btime TestParticle.solve($prob; dt=1e-3, trajectories=100)

param_td = prepare(E_field, B_field_td; species=Proton)
prob_td = TraceProblem(u0, tspan, param_td)

println("Benchmarking time-dependent field with current implementation:")
@btime TestParticle.solve($prob_td; dt=1e-3, trajectories=100)
