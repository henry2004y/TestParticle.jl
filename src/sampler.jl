# Particle sampling

using VelocityDistributionFunctions: Maxwellian, BiMaxwellian, Kappa, BiKappa, VelocityDistribution
import VelocityDistributionFunctions

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
     BiKappa(B::Vector{U}, u0::Vector{T}, ppar, pperp, n, kappa; m=mᵢ)

Construct a `BiKappa` distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, number density `n`, and
spectral index `kappa` in SI units. The default particle is proton.
"""
function VelocityDistributionFunctions.BiKappa(
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
   vpar = get_thermal_speed(ppar, n, m)
   vperp = get_thermal_speed(pperp, n, m)

   BiKappa(vperp, vpar, kappa, u0, B)
end

get_thermal_speed(p, n, m) = √(2 * p / (n * m))
