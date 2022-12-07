module TestParticle

using LinearAlgebra: norm, ×, ⋅
using Meshes: coordinates, spacing, embeddim, CartesianGrid
using Interpolations: interpolate, extrapolate, BSpline, Cubic, Line, OnGrid, Periodic,
   scale
using Distributions: Normal, MvNormal
using StaticArrays
using SnoopPrecompile

export prepare, sample
export trace!, trace_relativistic!, trace_normalized!, trace, trace_relativistic
export Proton, Electron, Ion, User
export Maxwellian, BiMaxwellian

include("utility/utility.jl")

"""
Type for the particles, `Proton`, `Electron`, `Ion`, or `User`.
"""
@enum Species Proton Electron Ion User

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

function getchargemass(species::Species, q::AbstractFloat, m::AbstractFloat)
   if species == Proton
      q = qᵢ
      m = mᵢ
   elseif species == Electron
      q = qₑ
      m = mₑ
   end
   q, m
end

function makegrid(grid::CartesianGrid{3, T}) where T
   gridmin = coordinates(minimum(grid))
   gridmax = coordinates(maximum(grid))
   Δx = spacing(grid)

   gridx = range(gridmin[1], gridmax[1], step=Δx[1])
   gridy = range(gridmin[2], gridmax[2], step=Δx[2])
   gridz = range(gridmin[3], gridmax[3], step=Δx[3])

   gridx, gridy, gridz
end

function makegrid(grid::CartesianGrid{2, T}) where T
   gridmin = coordinates(minimum(grid))
   gridmax = coordinates(maximum(grid))
   Δx = spacing(grid)

   gridx = range(gridmin[1], gridmax[1], step=Δx[1])
   gridy = range(gridmin[2], gridmax[2], step=Δx[2])

   gridx, gridy
end

function getinterp(A, gridx, gridy, gridz)
   @assert size(A,1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"

   Ax = @view A[1,:,:,:]
   Ay = @view A[2,:,:,:]
   Az = @view A[3,:,:,:]

   itp = extrapolate(interpolate(Ax, BSpline(Cubic(Line(OnGrid())))), NaN)
   interpx = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Ay, BSpline(Cubic(Line(OnGrid())))), NaN)
   interpy = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Az, BSpline(Cubic(Line(OnGrid())))), NaN)
   interpz = scale(itp, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:3]

      return SA[interpx(r...), interpy(r...), interpz(r...)]
   end

   return Field(get_field)
end

function getinterp(A, gridx, gridy)
   @assert size(A,1) == 3 && ndims(A) == 3 "Inconsistent 2D force field and grid!"

   Ax = @view A[1,:,:]
   Ay = @view A[2,:,:]
   Az = @view A[3,:,:]

   # The most common boundary condition for 2D is periodic.
   itp = extrapolate(interpolate(Ax, BSpline(Cubic(Periodic(OnGrid())))), Periodic())
   interpx = scale(itp, gridx, gridy)

   itp = extrapolate(interpolate(Ay, BSpline(Cubic(Periodic(OnGrid())))), Periodic())
   interpy = scale(itp, gridx, gridy)

   itp = extrapolate(interpolate(Az, BSpline(Cubic(Periodic(OnGrid())))), Periodic())
   interpz = scale(itp, gridx, gridy)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:2]

      return SA[interpx(r...), interpy(r...), interpz(r...)]
   end

   return Field(get_field)
end

# judge whether the field function is time dependent
function is_time_dependent(f)
   itp_param_number = 2
   numargs = [length(m.sig.parameters)-1 for m in methods(f)]
   itp = any(x -> x == itp_param_number, numargs)
   if !itp
      if any(x -> x == itp_param_number-1, numargs)
         return false
      else
         name = nameof(f)
         throw(ArgumentError(
            "All methods for the field function $name had too many arguments."))
      end
   else
      return true
   end
end

abstract type AbstractField{itd} end

"""
    Field{itd, F} <: AbstractField{itd}

A representation of a field function `f`, defined by:

time-independent field
```math
\\mathbf{F} = F(\\mathbf{x}),
```

time-dependent field
```math
\\mathbf{F} = F(\\mathbf{x}, t).
```

# Arguments
- `field_function::Function`: the function of field.
- `itd::Bool`: whether the field function is time dependent.
- `F`: the type of `field_function`.
"""
struct Field{itd, F} <: AbstractField{itd}
   field_function::F
   Field{itd, F}(field_function::F) where {itd, F} =
      isa(itd, Bool) ? new(field_function) : throw(ArgumentError("itd must be a boolean."))
end

Field(f::Function) = Field{is_time_dependent(f), typeof(f)}(f)

(f::AbstractField{true})(xu, t) = f.field_function(xu, t)
(f::AbstractField{true})(xu) = throw(ArgumentError("Time-dependent field function must have a time argument."))
(f::AbstractField{false})(xu, t) = f.field_function(xu)
(f::AbstractField{false})(xu) = f.field_function(xu)

"The type of parameter tuple for full test particle problem."
FullTPTuple = Tuple{Float64, Float64, AbstractField, AbstractField, AbstractField}

"The type of parameter tuple for normal test particle problem."
TPTuple = Tuple{Float64, Float64, AbstractField, AbstractField}

"The type of parameter tuple for normalized test particle problem."
TPNormalizedTuple = Tuple{AbstractFloat, AbstractField, AbstractField}

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

"""
    prepare(grid::CartesianGrid, E, B; species=Proton) -> (q, m, E, B)

Return a tuple consists of particle charge, mass for a prescribed `species` and interpolated
EM field functions.

    prepare(grid::CartesianGrid, E, B, B₀::Real; species=Proton) -> (Ω, E, B)

Return a tuple consists of gyrofrequency for normalization and interpolated EM field
functions given magnetic field scale B₀ in Tesla.

    prepare(grid::CartesianGrid, E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, interpolated EM field functions, and external force `F`.

    prepare(x::AbstractRange, y::AbstractRange, z::AbstractRange, E, B) -> (q, m, E, B)

Direct range input for uniform grid in 3D is also accepted.

    prepare(E, B; species=Proton, q=1.0, m=1.0) -> (q, m, E, B)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m` and analytic EM field functions. Prescribed `species` are `Electron` and
`Proton`; other species can be manually specified with `species=Ion/User`, `q` and `m`.

    prepare(E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, analytic EM field functions, and external force `F`.
"""
function prepare(grid::CartesianGrid, E::TE, B::TB; species::Species=Proton,
   q::AbstractFloat=1.0, m::AbstractFloat=1.0) where {TE, TB}

   q, m = getchargemass(species, q, m)

   if embeddim(grid) == 3
      gridx, gridy, gridz = makegrid(grid)
      E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz) : Field(E)
      B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz) : Field(B)
   elseif embeddim(grid) == 2
      gridx, gridy = makegrid(grid)
      E = TE <: AbstractArray ? getinterp(E, gridx, gridy) : Field(E)
      B = TB <: AbstractArray ? getinterp(B, gridx, gridy) : Field(B)
   end

   q, m, E, B
end

function prepare(grid::CartesianGrid, E::TE, B::TB, B₀::Real; species::Species=Proton,
   q::AbstractFloat=1.0, m::AbstractFloat=1.0) where {TE, TB}

   q, m = getchargemass(species, q, m)
   Ω = q*B₀/m

   if embeddim(grid) == 3
      gridx, gridy, gridz = makegrid(grid)
      E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz) : Field(E)
      B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz) : Field(B)
   elseif embeddim(grid) == 2
      gridx, gridy = makegrid(grid)
      E = TE <: AbstractArray ? getinterp(E, gridx, gridy) : Field(E)
      B = TB <: AbstractArray ? getinterp(B, gridx, gridy) : Field(B)
   end

   Ω, E, B
end

function prepare(grid::CartesianGrid, E::TE, B::TB, F::TF; species::Species=Proton,
   q::AbstractFloat=1.0, m::AbstractFloat=1.0) where {TE, TB, TF}

   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)

   E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz) : Field(E)
   B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz) : Field(B)
   F = TF <: AbstractArray ? getinterp(F, gridx, gridy, gridz) : Field(F)

   q, m, E, B, F
end

function prepare(x::AbstractRange, y::AbstractRange, z::AbstractRange, E::TE, B::TB;
   species::Species=Proton, q::AbstractFloat=1.0, m::AbstractFloat=1.0) where {TE, TB}

   q, m = getchargemass(species, q, m)

   E = TE <: AbstractArray ? getinterp(E, x, y, z) : Field(E)
   B = TB <: AbstractArray ? getinterp(B, x, y, z) : Field(B)

   q, m, E, B
end

function prepare(E, B; species::Species=Proton, q::AbstractFloat=1.0, m::AbstractFloat=1.0)
   q, m = getchargemass(species, q, m)
   B = Field(B)
   E = Field(E)

   q, m, E, B
end

function prepare(E, B, F; species::Species=Proton, q::AbstractFloat=1.0,
   m::AbstractFloat=1.0)

   q, m = getchargemass(species, q, m)
   B = Field(B)
   E = Field(E)
   F = Field(F)

   q, m, E, B, F
end

"""
    trace!(dy, y, p::TPTuple, t)
    trace!(dy, y, p::FullTPTuple, t)

ODE equations for charged particle moving in static EM field with in-place form.

ODE equations for charged particle moving in static EM field and external force field with
in-place form.
"""
function trace!(dy, y, p::TPTuple, t)
   q, m, E, B = p
   v = @view y[4:6]
   dy[1:3] = v
   dy[4:6] = q/m*(E(y, t) + v × (B(y, t)))
end

function trace!(dy, y, p::FullTPTuple, t)
   q, m, E, B, F = p
   v = @view y[4:6]
   dy[1:3] = v
   dy[4:6] = (q*(E(y, t) + v × (B(y, t))) + F(y, t)) / m
end

"""
    trace(y, p::TPTuple, t) -> SVector{6, Float64}
    trace(y, p::FullTPTuple, t) -> SVector{6, Float64}

ODE equations for charged particle moving in static EM field with out-of-place form.

ODE equations for charged particle moving in static EM field and external force field with
out-of-place form.
"""
function trace(y, p::TPTuple, t)
   q, m, E, B = p
   v = @view y[4:6]
   dx, dy, dz = v
   dux, duy, duz = q/m*(E(y, t) + v × (B(y, t)))
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

function trace(y, p::FullTPTuple, t)
   q, m, E, B, F = p
   v = @view y[4:6]
   dx, dy, dz = v
   dux, duy, duz = (q*(E(y, t) + v × (B(y, t))) + F(y, t)) / m
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

const FTLError = """
Particle faster than the speed of light!

If the initial velocity is slower than light and
adaptive timestepping of the solver is turned on, it
is better to set a small initial stepsize (dt) or
maximum dt for adaptive timestepping (dtmax).

More details about the keywords of initial stepsize
can be found in this documentation page:
https://diffeq.sciml.ai/stable/basics/common_solver_opts/#Stepsize-Control
"""

"""
    trace_relativistic!(dy, y, p::TPTuple, t)

ODE equations for relativistic charged particle moving in static EM field with in-place
form.
"""
function trace_relativistic!(dy, y, p::TPTuple, t)
   q, m, E, B = p

   u2 = y[4]^2 + y[5]^2 + y[6]^2
   c2 = c^2
   if u2 ≥ c2
      throw(DomainError(u2, FTLError))
   end

   γInv = √(1.0 - u2/c2)
   v = @view y[4:6]
   dy[1:3] = v
   dy[4:6] = q/m*γInv^3*(E(y, t) + v × (B(y, t)))
end

"""
    trace_relativistic(y, p::TPTuple, t) -> SVector{6, Float64}

ODE equations for relativistic charged particle moving in static EM field with out-of-place
form.
"""
function trace_relativistic(y, p::TPTuple, t)
   q, m, E, B = p

   u2 = y[4]^2 + y[5]^2 + y[6]^2
   c2 = c^2
   if u2 ≥ c2
      throw(DomainError(u2, FTLError))
   end

   γInv = √(1.0 - u2/c2)
   v = @view y[4:6]
   dx, dy, dz = v
   dux, duy, duz = q/m*γInv^3*(E(y, t) + v × (B(y, t)))
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

"""
    trace_normalized!(dy, y, p::TPNormalizedTuple, t)

Normalized ODE equations for charged particle moving in static EM field with in-place form.
If the field is in 2D X-Y plane, periodic boundary should be applied for the field in z via
the extrapolation function provided by Interpolations.jl.
"""
function trace_normalized!(dy, y, p::TPNormalizedTuple, t)
   Ω, E, B = p
   v = @view y[4:6]

   dy[1:3] = v
   dy[4:6] = Ω*(E(y, t) + v × B(y, t))
end

"""
    guiding_center(xu, param::Union{TPTuple, FullTPTuple})

Calculate the coordinates of the guiding center according to the phase space coordinates of a particle.
Reference: https://en.wikipedia.org/wiki/Guiding_center

A simple definition:
```math
\\mathbf{X}=\\mathbf{x}-m\\frac{\\mathbf{v}\\times\\mathbf{B}}{qB}
```
"""
function guiding_center(xu, param::Union{TPTuple, FullTPTuple})
   q, m, _, B_field = param
   t = xu[end]
   v = @view xu[4:6]
   Bv = B_field(xu, t)
   B = hypot(Bv...)
   # unit vector along B
   b = Bv./B
   # the vector of Larmor radius
   ρ = (b×v)./(q*B/m)
   X = xu[1:3] - ρ
   return X
end

"""
    get_gc(param::Union{TPTuple, FullTPTuple})

Get three function for plotting the orbit of guiding center.

For example:
```julia
param = prepare(E, B; species=Proton)
gc = get_gc(params)
# The definitions of stateinit, tspan, E and B are ignored.
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern7(); dt=2e-11)

orbit(sol, vars=gc)
```
"""
function get_gc(param::Union{TPTuple, FullTPTuple})
   gc_x(xu) = getindex(guiding_center(xu, param), 1)
   gc_y(xu) = getindex(guiding_center(xu, param), 2)
   gc_z(xu) = getindex(guiding_center(xu, param), 3)
   # Because of the design of the keyword 'vars', three coordinates must be divided into three functions.
   return (gc_x, gc_y, gc_z)
end

include("precompile.jl")

end