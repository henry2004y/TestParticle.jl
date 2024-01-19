# ---
# title: Ensemble tracing with extra saving
# id: demo_output_func
# date: 2023-12-20
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.4
# description: Tracing multiple charged particles in a static EM field
# ---

# This example demonstrates tracing multiple protons in an analytic E field and numerical B field.
# See [Demo: single tracing with additional diagnostics](@ref demo_savingcallback) for explaining the unit conversion.
# Also check [demo_ensemble](@ref demo_ensemble) for basic usages of the ensemble problem.
# The `output_func` argument can be used to change saving outputs. It works as a reduction function, but here we demonstrate how to add additional outputs.
# Besides the regular outputs, we also save the magnetic field along the trajectory, together with the parallel velocity.

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using Statistics
using LinearAlgebra
using Random
using CairoMakie
CairoMakie.activate!(type = "png")

seed = 1 # seed for random number
Random.seed!(seed)

"Set initial state for EnsembleProblem."
function prob_func(prob, i, repeat)
   B0 = prob.p[3](prob.u0)
   B0 = normalize(B0)

   Bperp1 = SA[0.0, -B0[3], B0[2]] |> normalize
   Bperp2 = B0 × Bperp1 |> normalize

   ## initial azimuthal angle
   ϕ = 2π*rand()
   ## initial pitch angle
   θ = acos(0.5)

   sinϕ, cosϕ = sincos(ϕ)
   prob.u0[4:6] = @. (B0*cos(θ) + Bperp1*(sin(θ)*cosϕ) + Bperp2*(sin(θ)*sinϕ)) * U₀

   prob
end

## Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 8
## Spatial extent along each dimension
x = range(-0.5, 0.5, length=nx)
y = range(-0.5, 0.5, length=ny)
z = range(-0.5, 0.5, length=nz)

## Numerical magnetic field
B = Array{Float32, 4}(undef, 3, nx, ny, nz)

B[1,:,:,:] .= 0.0
B[2,:,:,:] .= 0.0
B[3,:,:,:] .= 2.0

## Reference values for unit conversions
const B₀ = let Bmag = @views hypot.(B[1,:,:,:], B[2,:,:,:], B[3,:,:,:])
   sqrt(mean(vec(Bmag) .^ 2))
end

const Ω = abs(qᵢ) * B₀ / mᵢ
const t₀ = 1 / Ω  # [s]
const U₀ = 1.0    # [m/s]
const l₀ = U₀ * t₀ # [m]
const E₀ = U₀*B₀ # [V/m]
## Factor to scale the spatial coordinates
lscale = 2.0
## Scale the coordinates to control the number of discrete B values encountered
## along a given trajectory. In this case B is uniform, so it won't affect the result.
x /= lscale
y /= lscale
z /= lscale

## For full EM problems, the normalization of E and B should be done separately.
B ./= B₀
E(x) = SA[0.0/E₀, 0.0/E₀, 0.0/E₀]

## By default User type assumes q=1, m=1
## bc=2 uses periodic boundary conditions
param = prepare(x, y, z, E, B; species=User, bc=2)

x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [1.0, 0.0, 0.0] # initial velocity [v₀]
stateinit = [x0..., u0...]
tspan = (0.0, 2π) # one averaged gyroperiod based on B₀

saveat = tspan[2] / 40 # save interval

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)

"Set customized outputs for the ensemble problem."
function output_func(sol, i)
   Bfunc = sol.prob.p[3]
   b = Bfunc.(sol.u)

   μ = [@views b[i] ⋅ sol[4:6, i] / hypot(b[i]...) for i in eachindex(sol)]

   (sol.u, b, μ), false
end

## Number of trajectories
trajectories = 2

ensemble_prob = EnsembleProblem(prob; prob_func, output_func, safetycopy=false)
sols = solve(ensemble_prob, Vern9(), EnsembleThreads(); trajectories, saveat)

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Proton trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data,
)

for i in eachindex(sols)
   xp = [s[1] for s in sols[i][1]]
   yp = [s[2] for s in sols[i][1]]
   zp = [s[3] for s in sols[i][1]]
   lines!(ax, xp, yp, zp, label="$i")
end

f = DisplayAs.PNG(f) #hide