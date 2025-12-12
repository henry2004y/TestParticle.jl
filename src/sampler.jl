# Particle sampling

using VelocityDistributionFunctions: Maxwellian, BiMaxwellian, Kappa, VelocityDistribution
import VelocityDistributionFunctions
import Distributions
using Distributions: Gamma, rand
using Random: AbstractRNG, default_rng

"""
     Maxwellian(u0::AbstractVector{T}, p, n; m=mᵢ)

Construct a `Maxwellian` distribution with bulk velocity `u0`, thermal pressure `p`, and
number density `n` in SI units. The default particle is proton.
"""
function VelocityDistributionFunctions.Maxwellian(u0::AbstractVector{T}, p, n; m = mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = get_thermal_speed(p, n, m)

   Maxwellian(vth, u0)
end

"""
     BiMaxwellian(B::Vector{U}, u0::Vector{T}, ppar, pperp, n; m=mᵢ)

Construct a `BiMaxwellian` distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, and number density `n` in
SI units. The default particle is proton.
"""
function VelocityDistributionFunctions.BiMaxwellian(
      B::AbstractVector{U},
      u0::AbstractVector{T},
      ppar,
      pperp,
      n;
      m = mᵢ
) where
      {T <: AbstractFloat, U <: AbstractFloat}
   @assert length(u0) == 3 && length(B) == 3 "The field vector must have length 3!"
   vpar = get_thermal_speed(ppar, n, m)
   vperp = get_thermal_speed(pperp, n, m)

   BiMaxwellian(vperp, vpar, u0, B)
end

"""
     Kappa(u0::AbstractVector{T}, p, n, kappa; m=mᵢ)

Construct a `Kappa` distribution with bulk velocity `u0`, thermal pressure `p`, number density
`n`, and spectral index `kappa` in SI units. The default particle is proton.
"""
function VelocityDistributionFunctions.Kappa(u0::AbstractVector{T}, p, n, kappa; m = mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = get_thermal_speed(p, n, m)

   Kappa(vth, kappa, u0)
end

"""
Type for BiKappa velocity distributions with respect to the magnetic field.
"""
struct BiKappa{V <: AbstractVector, T <: AbstractFloat, U <: AbstractFloat} <: VelocityDistribution{T}
   "Unit magnetic field"
   b0::V
   "Bulk velocity"
   u0::V
   "Parallel thermal velocity"
   vth_para::T
   "Perpendicular thermal velocity"
   vth_perp::T
   "Spectral index"
   κ::T
end

"""
     BiKappa(B::Vector{U}, u0::Vector{T}, ppar, pperp, n, kappa; m=mᵢ)

Construct a `BiKappa` distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, number density `n`, and
spectral index `kappa` in SI units. The default particle is proton.
"""
function BiKappa(
      B::AbstractVector{U},
      u0::AbstractVector{T},
      ppar,
      pperp,
      n,
      kappa;
      m = mᵢ
) where
      {T <: AbstractFloat, U <: AbstractFloat}
   @assert length(u0) == 3 && length(B) == 3 "The field vector must have length 3!"
   b0 = normalize(B)
   vpar = get_thermal_speed(ppar, n, m)
   vperp = get_thermal_speed(pperp, n, m)

   BiKappa{typeof(B), T, U}(b0, u0, vpar, vperp, kappa)
end

"""
    rand(rng::AbstractRNG, d::BiKappa)

Sample a 3D velocity from a [`BiKappa`](@ref) distribution `d`.
"""
function Distributions._rand!(rng::AbstractRNG, vdf::BiKappa{V, T, U}, v::AbstractVector{<:Real}) where {V, T, U}
   # Sample Y ~ N(0, I)
   r1, r2, θ, ϕ = rand(rng, SVector{4, T})
   m1 = √(-2*log(r1))
   m2 = √(-2*log(r2))

   ypar = m1 * cospi(2θ)
   yperp1, yperp2 = m2 .* sincospi(2ϕ)

   # Sample U ~ Gamma(κ - 0.5, 2)
   u = rand(rng, Gamma(vdf.κ - 0.5, 2.0))

   # Scale factor
   scale_factor = sqrt(2 * (vdf.κ - 1.5) / u)

   vpar = ypar * scale_factor * vdf.vth_para
   vperp1 = yperp1 * scale_factor * vdf.vth_perp
   vperp2 = yperp2 * scale_factor * vdf.vth_perp

   tmp = vdf.b0 × SVector(1, 0, 0)
   bperp1 = if (tmp[1]^2 + tmp[2]^2 + tmp[3]^2) < 1e-6
      vdf.b0 × SVector(0, 1, 0)
   else
      tmp
   end
   bperp2 = vdf.b0 × bperp1

   v .= @. vdf.u0 + vpar * vdf.b0 + vperp1 * bperp1 + vperp2 * bperp2
end

function Base.show(io::IO, vdf::BiKappa)
   println(io, "BiKappa distribution")
   println(io, "B field direction:", vdf.b0)
   println(io, "Bulk velocity: ", vdf.u0)
   println(io, "Parallel thermal speed: ", vdf.vth_para)
   println(io, "Perpendicular thermal speed: ", vdf.vth_perp)
   println(io, "Spectral index: ", vdf.κ)
end

get_thermal_speed(p, n, m) = √(2 * p / (n * m))
