# Tracing multiple charged particles in a static EM field.
# This example takes advantage of the multithreading support in the ODE solver.
# A multiproc version is also available. Check the official documentation of
# DifferentialEquations.jl for details.
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using Meshes
using OrdinaryDiffEq
using Plots

"Initial state perturbation for EnsembleProblem."
function prob_func(prob, i, repeat)
   remake(prob, u0=rand()*prob.u0)
end

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

mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
   (x[1], y[1], z[1]),
   (Δx, Δy, Δz))

## Initialize particles

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]
   
param = prepare(mesh, E, B, species=Electron)
tspan = (0.0, 15.0)

trajectories = 10

## Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sol = solve(ensemble_prob, Tsit5(), EnsembleThreads();
   trajectories=trajectories, save_idxs=[1,2,3])

## Visualization

xyz = plot(sol, vars=(1,2,3), xlabel="x", ylabel="y", zlabel="z",
   aspect_ratio=:equal, label="electron", lw=2)