# Particle sampling

import VelocityDistributionFunctions

# Wrapper functions to avoid type piracy when adding convenience constructors

"""
     Maxwellian(args...; kw...)

Construct a `Maxwellian` distribution. Forwards to `VelocityDistributionFunctions.Maxwellian`.

     Maxwellian(u0::AbstractVector{T}, p, n; m=mᵢ)

Construct a `Maxwellian` distribution with bulk velocity `u0`, thermal pressure `p`, and
number density `n` in SI units. The default particle is proton.
"""
Maxwellian(args...; kwargs...) = VelocityDistributionFunctions.Maxwellian(args...; kwargs...)

function Maxwellian(u0::AbstractVector{T}, p, n; m = mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = get_thermal_speed(p, n, m)

   VelocityDistributionFunctions.Maxwellian(vth; u0)
end

"""
     BiMaxwellian(args...; kw...)

Construct a `BiMaxwellian` distribution. Forwards to `VelocityDistributionFunctions.BiMaxwellian`.

     BiMaxwellian(B::Vector{U}, u0::Vector{T}, ppar, pperp, n; m=mᵢ)

Construct a `BiMaxwellian` distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, and number density `n` in
SI units. The default particle is proton.
"""
BiMaxwellian(args...; kwargs...) = VelocityDistributionFunctions.BiMaxwellian(args...; kwargs...)

function BiMaxwellian(
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

   VelocityDistributionFunctions.BiMaxwellian(vperp, vpar, B; u0)
end

"""
     Kappa(args...; kw...)

Construct a `Kappa` distribution. Forwards to `VelocityDistributionFunctions.Kappa`.

     Kappa(u0::AbstractVector{T}, p, n, kappa; m=mᵢ)

Construct a `Kappa` distribution with bulk velocity `u0`, thermal pressure `p`, number density
`n`, and spectral index `kappa` in SI units. The default particle is proton.
"""
Kappa(args...; kwargs...) = VelocityDistributionFunctions.Kappa(args...; kwargs...)

function Kappa(u0::AbstractVector{T}, p, n, kappa; m = mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = get_thermal_speed(p, n, m)

   VelocityDistributionFunctions.Kappa(vth, kappa; u0)
end

"""
     BiKappa(args...; kw...)

Construct a `BiKappa` distribution. Forwards to `VelocityDistributionFunctions.BiKappa`.

     BiKappa(B::Vector{U}, u0::Vector{T}, ppar, pperp, n, kappa; m=mᵢ)

Construct a `BiKappa` distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, number density `n`, and
spectral index `kappa` in SI units. The default particle is proton.
"""
BiKappa(args...; kwargs...) = VelocityDistributionFunctions.BiKappa(args...; kwargs...)

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
   vpar = get_thermal_speed(ppar, n, m)
   vperp = get_thermal_speed(pperp, n, m)

   VelocityDistributionFunctions.BiKappa(vperp, vpar, kappa, B; u0)
end

get_thermal_speed(p, n, m) = √(2 * p / (n * m))
