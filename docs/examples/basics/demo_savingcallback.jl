# ---
# title: Single tracing with additional diagnostics
# id: demo_savingcallback
# date: 2023-12-17
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.4
# description: Tracing one charged particle with additional diagnostics
# ---

# This example demonstrates tracing one proton in an analytic E field and numerical B field.
# It also combines one type of normalization using a reference velocity `U₀`, a reference magnetic field `B₀`, and a reference time `1/Ω`, where `Ω` is the gyrofrequency.
# This indicates that in the dimensionless units, a proton with initial perpendicular velocity 1 will possess a gyro-radius of 1.
# In the dimensionless spatial coordinates, we can zoom in/out the EM field to control the number of discrete points encountered in a gyroperiod.
# For example, if `dx=dy=dz=1`, it means that a particle with perpendicular velocity 1 will "see" one discrete point along a certain direction oriented from the gyro-center within the gyro-radius;
# if `dx=dy=dz=0.5`, then the particle will "see" two discrete points.
# MHD models, for instance, are dimensionless by nature.
# If we simulate a turbulence field with MHD, we want to include more discrete points within a gyro-radius for the effect of small scale perturbations to take place. (Otherwise within one gyro-period all you will see is a nice-looking helix!)
# However, we cannot simply shrink the spatial coordinates as we wish, otherwise we will quickly encounter the boundary of our simulation.

# The `SavingCallback` from DiffEqCallbacks.jl can be used to save additional outputs for diagnosis. Here we save the magnetic field along the trajectory, together with the parallel velocity.
# Note that `SavingCallback` is currently not compatible with ensemble problems; for multiple particle tracing with customized outputs, see [Demo: ensemble tracing with extra saving](@ref demo_output_func).

using TestParticle
using TestParticle: qᵢ, mᵢ
using OrdinaryDiffEq
using StaticArrays
using Statistics
using LinearAlgebra
using DiffEqCallbacks

## Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 8
## Numerical magnetic field given in customized dimensionless units
B = Array{Float32, 4}(undef, 3, nx, ny, nz)

B[1,:,:,:] .= 0.0
B[2,:,:,:] .= 0.0
B[3,:,:,:] .= 2.0

## Reference values for unit conversions between the customized and default dimensionless
## units
const B₀ = 2.0 
const Ω = B₀
const t₀ = 1 / Ω
const U₀ = 1.0
const l₀ = U₀ * t₀
const E₀ = U₀ * B₀

## Convert from customized to default dimensionless units
## Dimensionless spatial extent along each dimension with resolution 1
x = range(0, nx-1, length=nx)
y = range(0, ny-1, length=ny)
z = range(0, nz-1, length=nz)
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

tspan = (0.0, π) # half averaged gyroperiod based on B₀

## Dummy initial state; positions have units l₀; velocities have units U₀
stateinit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)

prob.u0[4:6] = let
   B0 = prob.p[3](prob.u0)
   B0 = normalize(B0)

   Bperp1 = SA[0.0, -B0[3], B0[2]] |> normalize
   Bperp2 = B0 × Bperp1 |> normalize

   ## initial velocity azimuthal angle
   ϕ = 2π*0
   ## initial velocity pitch angle w.r.t. B
   θ = acos(0.0)

   sinϕ, cosϕ = sincos(ϕ)
   @. (B0*cos(θ) + Bperp1*(sin(θ)*cosϕ) + Bperp2*(sin(θ)*sinϕ)) * U₀
end

saved_values = SavedValues(Float64, Tuple{SVector{3, Float64}, Float64})

function save_B_mu(u, t, integrator)
   b = integrator.p[3](u)
   μ = @views b ⋅ u[4:6] / hypot(b...)

   b, μ
end

cb = SavingCallback(save_B_mu, saved_values)

sol = solve(prob, Vern9(); callback=cb)

# The extra values are saved in `saved_values`:

saved_values