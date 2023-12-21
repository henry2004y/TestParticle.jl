# ---
# title: Single tracing with additional diagnostics
# id: demo_savingcallback
# date: 2023-12-17
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.4
# description: Tracing one charged particle with additional diagnostics
# ---

# This example demonstrates tracing one proton in an analytic E field and numerical B field.
# It also combines one type of normalization using a reference velocity `U₀`, a reference magnetic field `B₀`, and an initial reference gyroradius `rL`.
# One way to think of this is that a proton with initial perpendicular velocity `U₀` will gyrate with a radius `rL`.
# The `SavingCallback` from DiffEqCallbacks.jl can be used to save additional outputs for diagnosis. Here we save the magnetic field along the trajectory, together with the parallel velocity.
# Note that `SavingCallback` is currently not compatible with ensemble problems; for multiple particle tracing with customized outputs, see [demo_output_func](@ref demo_output_func).

using TestParticle
using OrdinaryDiffEq
using StaticArrays
using Statistics
using LinearAlgebra
using DiffEqCallbacks

## Analytical electric field
E(x) = SA[0.0, 0.0, 0.0]
## Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 8
## Spatial extent along each dimension
x = range(-10.0, 10.0, length=nx)
y = range(-10.0, 10.0, length=ny)
z = range(-10.0, 10.0, length=nz)

## Numerical magnetic field
B = Array{Float32, 4}(undef, 3, nx, ny, nz)

B[1,:,:,:] .= 0.0
B[2,:,:,:] .= 0.0
B[3,:,:,:] .= 1.0

Bmag = @views hypot.(B[1,:,:,:], B[2,:,:,:], B[3,:,:,:])

B₀ = sqrt(mean(vec(Bmag) .^ 2))

const U₀ = 1.0
const rL = 4.0

Ω = q2mc = U₀ / (rL*B₀)

## By default User type assumes q=1, m=1
## bc=2 uses periodic boundary conditions
param = prepare(x, y, z, E, B, Ω; species=User, bc=2)

tspan = (0.0, 2π/Ω) # one averaged gyroperiod based on B₀

## Dummy initial state
stateinit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)

prob.u0[4:6] = let
   x0 = [0.0, 0.0, 0.0] # initial position [l₀]

   B0 = prob.p[3](prob.u0)
   B0 = normalize(B0)

   Bperp1 = SA[0.0, -B0[3], B0[2]] |> normalize
   Bperp2 = B0 × Bperp1 |> normalize

   ## initial azimuthal angle
   ϕ = 2π*rand()
   ## initial pitch angle
   θ = acos(0.5)

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