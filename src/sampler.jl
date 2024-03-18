# Particle sampling

"""
Abstract type for velocity distribution functions.
"""
abstract type VDF end

"""
Type for Maxwellian velocity distributions. 
"""
struct Maxwellian{T<:AbstractFloat} <: VDF
   "Bulk velocity"
   u0::Vector{T}
   "Thermal speed"
   vth::T
end

"""
    Maxwellian(u0::Vector{T}, p::T, n; m=mᵢ)

Construct a Maxwellian distribution with bulk velocity `u0`, thermal pressure `p`, and
number density `n` in SI units. The default particle is proton.
"""
function Maxwellian(u0::Vector{T}, p::T, n; m=mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = √(p / (n * m))

   Maxwellian{T}(u0, vth)
end

"""
Type for BiMaxwellian velocity distributions with respect to the magnetic field.
"""
struct BiMaxwellian{T<:AbstractFloat, U} <: VDF
   "Unit magnetic field"
   b0::Vector{U}
   "Bulk velocity"
   u0::Vector{T}
   "Parallel thermal velocity"
   vthpar::T
   "Perpendicular thermal velocity"
   vthperp::T
end

"""
    BiMaxwellian(B::Vector{U}, u0::Vector{T}, ppar::T, pperp::T, n; m=mᵢ)

Construct a BiMaxwellian distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, and number density `n` in
SI units. The default particle is proton.
"""
function BiMaxwellian(B::Vector{U}, u0::Vector{T}, ppar::T, pperp::T, n; m=mᵢ) where
   {T <: AbstractFloat, U <: AbstractFloat}
   @assert length(u0) == 3 && length(B) == 3 "The field vector must have length 3!"
   b0 = normalize(B)
   vpar = √(ppar / (n * m))
   vperp = √(pperp / (n * m))

   BiMaxwellian{T, U}(b0, u0, vpar, vperp)
end

"""
    sample(vdf::Maxwellian)

Sample a 3D velocity from a [`Maxwellian`](@ref) distribution `vdf` using the Box-Muller method.

    sample(vdf::BiMaxwellian)

Sample a 3D velocity from a [`BiMaxwellian`](@ref) distribution `vdf` using the Box-Muller method.
"""
function sample(vdf::Maxwellian{T}) where T
   r1, r2, θ, ϕ = rand(T, 4)
   m1 = √(-2*log(r1))
   m2 = √(-2*log(r2))

   v1 = vdf.vth * m1 * cospi(2θ)
   v2, v3 = vdf.vth .* m2 .* sincospi(2ϕ)
   v1v = v1 .* SVector(1, 0, 0)
   v2v = v2 .* SVector(0, 1, 0)
   v3v = v3 .* SVector(0, 0, 1)

   v = @. vdf.u0 + v1v + v2v + v3v
end

function sample(vdf::BiMaxwellian{T, U}) where {T, U}
   r1, r2, θ, ϕ = rand(T, 4)
   m1 = √(-2*log(r1))
   m2 = √(-2*log(r2))

   vpar = vdf.vthpar * m1 * cospi(2θ)
   vperp1, vperp2 = vdf.vthperp .* m2 .* sincospi(2ϕ)

   tmp = vdf.b0 × SVector(1, 0, 0)
   bperp1 =
      if hypot(tmp...) < 1e-3
         vdf.b0 × SVector(0, 1, 0)
      else
         tmp
      end
   bperp2 = vdf.b0 × bperp1

   v = @. vdf.u0 + vpar * vdf.b0  + vperp1 * bperp1 + vperp2 * bperp2
end

function Base.show(io::IO, vdf::Maxwellian)
   println(io, "Isotropic Maxwellian distribution")
   println(io, "Bulk velocity: ", vdf.u0)
   println(io, "Thermal speed: ", vdf.vth)
end

function Base.show(io::IO, vdf::BiMaxwellian)
   println(io, "BiMaxwellian distribution")
   println(io, "B field direction:", vdf.b0)
   println(io, "Bulk velocity: ", vdf.u0)
   println(io, "Parallel thermal speed: ", vdf.vthpar)
   println(io, "Perpendicular thermal speed: ", vdf.vthperp)
end