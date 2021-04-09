# Tracing charged particle in a static EM field.
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using Meshes, DifferentialEquations
using Plots

## Initialize grid and field

x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
z = range(-10, 10, length=25)
B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

B[3,:,:,:] .= 1e-11
E[3,:,:,:] .= 5e-13

Δx = x[2] - x[1]
Δy = y[2] - y[1]
Δz = z[2] - z[1]

grid = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
   (x[1], y[1], z[1]),
   (Δx, Δy, Δz))

## Initialize particles

isAnalytic = false
trajectories = 1

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param_electron = prepare(grid, E, B, species="electron")
tspan_electron = (0.0, 15.0)

param_proton = prepare(grid, E, B, species="proton")
tspan_proton = (0.0, 10.0)

## Solve for the trajectories

prob_e = ODEProblem(trace_numeric!, stateinit, tspan_electron, param_electron)
prob_p = ODEProblem(trace_numeric!, stateinit, tspan_proton, param_proton)

sol_e = solve(prob_e; save_idxs=[1,2,3], alg_hints=[:nonstiff])
sol_p = solve(prob_p; save_idxs=[1,2,3], alg_hints=[:nonstiff])

## Visualization

xyz = plot(sol_e, vars=(1,2,3), xlabel="x", ylabel="y", zlabel="z",
   aspect_ratio=:equal, label="electron", lw=2)

plot!(sol_p, vars=(1,2,3), label="proton", lw=2, camera = (70,50))
