using TestParticle
using TestParticle: qᵢ, mᵢ
using OrdinaryDiffEq
using StaticArrays
using Statistics
using LinearAlgebra
using DiffEqCallbacks

# Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 8
# Spatial coordinates given in customized units
x = range(-0.5, 0.5, length=nx)
y = range(-0.5, 0.5, length=ny)
z = range(-0.5, 0.5, length=nz)
# Numerical magnetic field given in customized units
B = Array{Float32, 4}(undef, 3, nx, ny, nz)

B[1,:,:,:] .= 0.0
B[2,:,:,:] .= 0.0
B[3,:,:,:] .= 2.0

# Reference values for unit conversions between the customized and dimensionless units
const B₀ = let Bmag = @views hypot.(B[1,:,:,:], B[2,:,:,:], B[3,:,:,:])
   sqrt(mean(vec(Bmag) .^ 2))
end
const U₀ = 1.0
const l₀ = 4*nx
const t₀ = l₀ / U₀
const E₀ = U₀ * B₀

### Convert from customized to default dimensionless units
# Dimensionless spatial extents [l₀]
x /= l₀
y /= l₀
z /= l₀
# For full EM problems, the normalization of E and B should be done separately.
B ./= B₀
E(x) = SA[0.0/E₀, 0.0/E₀, 0.0/E₀]

# By default User type assumes q=1, m=1; bc=2 uses periodic boundary conditions
param = prepare(x, y, z, E, B; species=User, bc=2)

tspan = (0.0, π) # half averaged gyroperiod based on B₀

# Dummy initial state; positions have units l₀; velocities have units U₀
stateinit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)

prob.u0[4:6] = let
   B0 = prob.p[3](prob.u0)
   B0 = normalize(B0)

   Bperp1 = SA[0.0, -B0[3], B0[2]] |> normalize
   Bperp2 = B0 × Bperp1 |> normalize

   # initial velocity azimuthal angle
   ϕ = 2π*0
   # initial velocity pitch angle w.r.t. B
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

saved_values

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl