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
   uth::T

   function Maxwellian(u0::Vector{T}, uth::T) where T
      @assert length(u0) == 3 "Bulk velocity must have length 3!"
      new{T}(u0, uth)
   end
end

"""
Type for BiMaxwellian velocity distributions.
"""
struct BiMaxwellian{T<:AbstractFloat, U} <: VDF
   "Bulk velocity"
   u0::Vector{T}
   "Parallel thermal speed"
   uthpar::T
   "Perpendicular thermal speed"
   uthperp::T
   "Unit magnetic field"
   b0::Vector{U}

   function BiMaxwellian(u0::Vector{T}, upar::T, uperp::T, B::Vector{U}) where
      {T <: AbstractFloat, U <: AbstractFloat}
      @assert length(u0) == 3 && length(B) == 3 "The field vector must have length 3!"
      b0 = B ./ hypot(B...)
      new{T, U}(u0, upar, uperp, b0)
   end
end


"""
    sample(vdf::Maxwellian, nparticles::Int)

Sample velocities from a [`Maxwellian`](@ref) distribution `vdf` with `npoints`.

    sample(vdf::BiMaxwellian, nparticles::Int)

Sample velocities from a [`BiMaxwellian`](@ref) distribution `vdf` with `npoints`.
"""
function sample(vdf::Maxwellian, nparticles::Int)
   sqr2 = typeof(vdf.uth)(√2)
   # Convert from thermal speed to std
   σ = vdf.uth / sqr2
   v = σ .* randn(typeof(vdf.uth), 3, nparticles) .+ vdf.u0
end

function sample(vdf::BiMaxwellian{T, U}, nparticles::Int) where {T, U}
   sqr2 = T(√2)
   # Convert from thermal speed to std
   σpar = vdf.uthpar / sqr2
   σperp = vdf.uthperp / sqr2
   # Transform to Cartesian grid
   v = fill(vdf.u0, (3, nparticles))
   vrand = σpar .* randn(T, nparticles)
   vrand = reshape(vrand, 1, nparticles)
   vpar = repeat(vrand, outer=3)
   @inbounds for i in 1:3, ip in 1:nparticles
      vpar[i,ip] = vpar[i,ip]*vdf.b0[i] + vdf.u0[i]
   end
   # Sample vectors on a 2D plane
   μ = zeros(SVector{2,T})
   σ = SA[σperp 0; 0 σperp]
   d = MvNormal(μ, σ)
   vrand = rand(d, nparticles)

   vperp = zeros(T, (3, nparticles))
   # Rotate vectors to be perpendicular to b̂
   k = SVector{3, T}(0, 0, 1)
   axis = vdf.b0 × k::SVector{3, T}
   θ = acos(vdf.b0 ⋅ k)
   R = get_rotation_matrix(axis, θ)
   @inbounds for ip in 1:nparticles
      vperp[:,ip] = R * SA[vrand[1,ip], vrand[2,ip], 0]
   end

   v = vpar .+ vperp
end