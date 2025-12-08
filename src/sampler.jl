# Particle sampling

using SpecialFunctions: gamma
"""
Abstract type for velocity distribution functions.
"""
abstract type VDF end

"""
Type for Maxwellian velocity distributions.
"""
struct Maxwellian{V <: AbstractVector, T <: AbstractFloat} <: VDF
   "Bulk velocity"
   u0::V
   "Thermal speed"
   vth::T
end

"""
     Maxwellian(u0::AbstractVector{T}, p, n; m=mᵢ)

Construct a Maxwellian distribution with bulk velocity `u0`, thermal pressure `p`, and
number density `n` in SI units. The default particle is proton.
"""
function Maxwellian(u0::AbstractVector{T}, p, n; m = mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = √(p / (n * m))

   Maxwellian{typeof(u0), T}(u0, vth)
end

"""
Type for BiMaxwellian velocity distributions with respect to the magnetic field.
"""
struct BiMaxwellian{V <: AbstractVector, T <: AbstractFloat, U <: AbstractFloat} <: VDF
   "Unit magnetic field"
   b0::V
   "Bulk velocity"
   u0::V
   "Parallel thermal velocity"
   vthpar::T
   "Perpendicular thermal velocity"
   vthperp::T
end

"""
     BiMaxwellian(B::Vector{U}, u0::Vector{T}, ppar, pperp, n; m=mᵢ)

Construct a BiMaxwellian distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, and number density `n` in
SI units. The default particle is proton.
"""
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
   b0 = normalize(B)
   vpar = √(ppar / (n * m))
   vperp = √(pperp / (n * m))

   BiMaxwellian{typeof(B), T, U}(b0, u0, vpar, vperp)
end

"""
Type for Kappa velocity distributions.
"""
struct Kappa{V <: AbstractVector, T <: AbstractFloat} <: VDF
   "Bulk velocity"
   u0::V
   "Thermal speed"
   vth::T
   "Spectral index"
   kappa::T
end

"""
     Kappa(u0::AbstractVector{T}, p, n, kappa; m=mᵢ)

Construct a Kappa distribution with bulk velocity `u0`, thermal pressure `p`, number density
`n`, and spectral index `kappa` in SI units. The default particle is proton.
"""
function Kappa(u0::AbstractVector{T}, p, n, kappa; m = mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = √(p / (n * m))

   Kappa{typeof(u0), T}(u0, vth, kappa)
end

"""
Type for BiKappa velocity distributions with respect to the magnetic field.
"""
struct BiKappa{V <: AbstractVector, T <: AbstractFloat, U <: AbstractFloat} <: VDF
   "Unit magnetic field"
   b0::V
   "Bulk velocity"
   u0::V
   "Parallel thermal velocity"
   vthpar::T
   "Perpendicular thermal velocity"
   vthperp::T
   "Spectral index"
   kappa::T
end

"""
     BiKappa(B::Vector{U}, u0::Vector{T}, ppar, pperp, n, kappa; m=mᵢ)

Construct a BiKappa distribution with magnetic field `B`, bulk velocity `u0`, parallel
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
   vpar = √(ppar / (n * m))
   vperp = √(pperp / (n * m))

   BiKappa{typeof(B), T, U}(b0, u0, vpar, vperp, kappa)
end

"""
Type for SelfSimilar velocity distributions (isotropic).
"""
struct SelfSimilar{V <: AbstractVector, T <: AbstractFloat} <: VDF
   "Bulk velocity"
   u0::V
   "Thermal speed"
   vth::T
   "Self-similar exponent"
   s::T
end

"""
     SelfSimilar(u0::AbstractVector{T}, p, n, s; m=mᵢ)

Construct a SelfSimilar distribution with bulk velocity `u0`, thermal pressure `p`, number density
`n`, and self-similar exponent `s` in SI units. The default particle is proton.
"""
function SelfSimilar(u0::AbstractVector{T}, p, n, s; m = mᵢ) where T
   @assert length(u0) == 3 "Bulk velocity must have length 3!"
   vth = √(p / (n * m))

   SelfSimilar{typeof(u0), T}(u0, vth, s)
end

"""
Type for BiSelfSimilar velocity distributions with respect to the magnetic field.
"""
struct BiSelfSimilar{V <: AbstractVector, T <: AbstractFloat, U <: AbstractFloat} <: VDF
   "Unit magnetic field"
   b0::V
   "Bulk velocity"
   u0::V
   "Parallel thermal velocity"
   vthpar::T
   "Perpendicular thermal velocity"
   vthperp::T
   "Parallel self-similar exponent"
   p_exp::T
   "Perpendicular self-similar exponent"
   q_exp::T
end

"""
     BiSelfSimilar(B::Vector{U}, u0::Vector{T}, ppar, pperp, n, p_exp, q_exp; m=mᵢ)

Construct a BiSelfSimilar distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, number density `n`, parallel
exponent `p_exp` and perpendicular exponent `q_exp` in SI units. The default particle is proton.
"""
function BiSelfSimilar(
      B::AbstractVector{U},
      u0::AbstractVector{T},
      ppar,
      pperp,
      n,
      p_exp,
      q_exp;
      m = mᵢ
) where
      {T <: AbstractFloat, U <: AbstractFloat}
   @assert length(u0) == 3 && length(B) == 3 "The field vector must have length 3!"
   b0 = normalize(B)
   vpar = √(ppar / (n * m))
   vperp = √(pperp / (n * m))

   BiSelfSimilar{typeof(B), T, U}(b0, u0, vpar, vperp, p_exp, q_exp)
end

"""
     sample(vdf::Maxwellian)

Sample a 3D velocity from a [`Maxwellian`](@ref) distribution `vdf` using the Box-Muller method.

     sample(vdf::BiMaxwellian)

Sample a 3D velocity from a [`BiMaxwellian`](@ref) distribution `vdf` using the Box-Muller method.

     sample(vdf::Kappa)

Sample a 3D velocity from a [`Kappa`](@ref) distribution `vdf`.

     sample(vdf::BiKappa)

Sample a 3D velocity from a [`BiKappa`](@ref) distribution `vdf`.

     sample(vdf::SelfSimilar)

Sample a 3D velocity from a [`SelfSimilar`](@ref) distribution `vdf`.

     sample(vdf::BiSelfSimilar)

Sample a 3D velocity from a [`BiSelfSimilar`](@ref) distribution `vdf`.
"""
function sample(vdf::Maxwellian{U, T}) where {U, T}
   r1, r2, θ, ϕ = rand(SVector{4, T})
   m1 = √(-2*log(r1))
   m2 = √(-2*log(r2))

   v1 = vdf.vth * m1 * cospi(2θ)
   v2, v3 = vdf.vth .* m2 .* sincospi(2ϕ)
   v1v = v1 .* SVector(1, 0, 0)
   v2v = v2 .* SVector(0, 1, 0)
   v3v = v3 .* SVector(0, 0, 1)

   v = @. vdf.u0 + v1v + v2v + v3v
end

function sample(vdf::BiMaxwellian{V, T, U}) where {V, T, U}
   r1, r2, θ, ϕ = rand(SVector{4, T})
   m1 = √(-2*log(r1))
   m2 = √(-2*log(r2))

   vpar = vdf.vthpar * m1 * cospi(2θ)
   vperp1, vperp2 = vdf.vthperp .* m2 .* sincospi(2ϕ)

   tmp = vdf.b0 × SVector(1, 0, 0)
   bperp1 = if (tmp[1]^2 + tmp[2]^2 + tmp[3]^2) < 1e-6
      vdf.b0 × SVector(0, 1, 0)
   else
      tmp
   end
   bperp2 = vdf.b0 × bperp1

   v = @. vdf.u0 + vpar * vdf.b0 + vperp1 * bperp1 + vperp2 * bperp2
end

function sample(vdf::Kappa{U, T}) where {U, T}
   # Sample Y ~ N(0, I)
   y = randn(SVector{3, T})

   # Sample U ~ Gamma(κ - 0.5, 2)
   u = rand_gamma(vdf.kappa - 0.5; scale = 2.0)

   # Scale factor
   # Variance of each component should be vth^2
   # Currently variance is vth^2 / 2 without the sqrt(2) factor
   scale_factor = sqrt(2 * (vdf.kappa - 1.5) * vdf.vth^2 / u)

   v = vdf.u0 + y * scale_factor
end

function sample(vdf::BiKappa{V, T, U}) where {V, T, U}
   # Sample Y ~ N(0, I)
   r1, r2, θ, ϕ = rand(SVector{4, T})
   m1 = √(-2*log(r1))
   m2 = √(-2*log(r2))

   ypar = m1 * cospi(2θ)
   yperp1, yperp2 = m2 .* sincospi(2ϕ)

   # Sample U ~ Gamma(κ - 0.5, 2)
   u = rand_gamma(vdf.kappa - 0.5; scale = 2.0)

   # Scale factor
   scale_factor = sqrt(2 * (vdf.kappa - 1.5) / u)

   vpar = ypar * scale_factor * vdf.vthpar
   vperp1 = yperp1 * scale_factor * vdf.vthperp
   vperp2 = yperp2 * scale_factor * vdf.vthperp

   tmp = vdf.b0 × SVector(1, 0, 0)
   bperp1 = if (tmp[1]^2 + tmp[2]^2 + tmp[3]^2) < 1e-6
      vdf.b0 × SVector(0, 1, 0)
   else
      tmp
   end
   bperp2 = vdf.b0 × bperp1

   v = @. vdf.u0 + vpar * vdf.b0 + vperp1 * bperp1 + vperp2 * bperp2
end

function sample(vdf::SelfSimilar{U, T}) where {U, T}
   # Sample 3 independent Generalized Normal variables
   # Parallel (z) and perpendicular (x, y) components
   # Here isotropic, so all same exponent s
   # Using Cartesian product definition consistent with normalization in some contexts,
   # but usually "SelfSimilar" implies spherical symmetry for isotropic.
   # If isotropic spherical: f(v) ~ exp(-(v/vth)^s).
   # Then v ~ GenGamma?
   # If we assume Cartesian separability (as in BiSelfSimilar normalization):
   # Calculate scale α such that variance is vth^2
   # α = vth * sqrt(Γ(1/s) / Γ(3/s))
   scale_val = vdf.vth * sqrt(gamma(1 / vdf.s) / gamma(3 / vdf.s))

   vx = rand_gen_normal(vdf.s; scale = scale_val)
   vy = rand_gen_normal(vdf.s; scale = scale_val)
   vz = rand_gen_normal(vdf.s; scale = scale_val)

   v = vdf.u0 + SVector(vx, vy, vz)
end

function sample(vdf::BiSelfSimilar{V, T, U}) where {V, T, U}
   # Parallel component: Generalized Normal with exponent p
   # α_par = vthpar * sqrt(Γ(1/p) / Γ(3/p))
   scale_par = vdf.vthpar * sqrt(gamma(1 / vdf.p_exp) / gamma(3 / vdf.p_exp))
   vpar = rand_gen_normal(vdf.p_exp; scale = scale_par)

   # Perpendicular component:
   # Assuming cylindrical symmetry for physical reasons in magnetized plasma.
   # f(vperp) ~ exp(-(vperp/α_perp)^q) * vperp
   # We want <vperp^2> = 2 * vthperp^2
   # <vperp^2> = α^2 * Γ(4/q) / Γ(2/q)
   # α^2 = 2 * vthperp^2 * Γ(2/q) / Γ(4/q)
   # α = vthperp * sqrt(2 * Γ(2/q) / Γ(4/q))
   scale_perp = vdf.vthperp * sqrt(2 * gamma(2 / vdf.q_exp) / gamma(4 / vdf.q_exp))

   # vperp = α * U^(1/q) where U ~ Gamma(2/q, 1)
   u_perp = rand_gamma(2.0 / vdf.q_exp; scale = 1.0)
   vperp_mag = scale_perp * u_perp^(1.0 / vdf.q_exp)
   # Random phase
   phi = 2π * rand()
   vperp1 = vperp_mag * cos(phi)
   vperp2 = vperp_mag * sin(phi)

   tmp = vdf.b0 × SVector(1, 0, 0)
   bperp1 = if (tmp[1]^2 + tmp[2]^2 + tmp[3]^2) < 1e-6
      vdf.b0 × SVector(0, 1, 0)
   else
      tmp
   end
   bperp2 = vdf.b0 × bperp1

   v = @. vdf.u0 + vpar * vdf.b0 + vperp1 * bperp1 + vperp2 * bperp2
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

function Base.show(io::IO, vdf::Kappa)
   println(io, "Isotropic Kappa distribution")
   println(io, "Bulk velocity: ", vdf.u0)
   println(io, "Thermal speed: ", vdf.vth)
   println(io, "Spectral index: ", vdf.kappa)
end

function Base.show(io::IO, vdf::BiKappa)
   println(io, "BiKappa distribution")
   println(io, "B field direction:", vdf.b0)
   println(io, "Bulk velocity: ", vdf.u0)
   println(io, "Parallel thermal speed: ", vdf.vthpar)
   println(io, "Perpendicular thermal speed: ", vdf.vthperp)
   println(io, "Spectral index: ", vdf.kappa)
end

function Base.show(io::IO, vdf::SelfSimilar)
   println(io, "Isotropic Self-Similar distribution")
   println(io, "Bulk velocity: ", vdf.u0)
   println(io, "Thermal speed: ", vdf.vth)
   println(io, "Exponent: ", vdf.s)
end

function Base.show(io::IO, vdf::BiSelfSimilar)
   println(io, "BiSelfSimilar distribution")
   println(io, "B field direction:", vdf.b0)
   println(io, "Bulk velocity: ", vdf.u0)
   println(io, "Parallel thermal speed: ", vdf.vthpar)
   println(io, "Perpendicular thermal speed: ", vdf.vthperp)
   println(io, "Parallel exponent: ", vdf.p_exp)
   println(io, "Perpendicular exponent: ", vdf.q_exp)
end
