module TestParticle

using LinearAlgebra: norm, ×
using Meshes
using Interpolations
using StaticArrays

export prepare, trace!, trace_relativistic!
export trace, trace_relativistic
export Proton, Electron, Ion, User

include("utility/utility.jl")

"""
Type for the particles, `Proton`, `Electron`, `Ion`, or `User`.
"""
@enum Species Proton Electron Ion User

function getchargemass(species::Species, q, m)
   if species == Proton
      q = qᵢ
      m = mᵢ
   elseif species == Electron
      q = qₑ
      m = mₑ
   end
   q, m
end

function makegrid(grid)
   gridmin = coordinates(minimum(grid))
   gridmax = coordinates(maximum(grid))
   Δx = spacing(grid)

   gridx = range(gridmin[1], gridmax[1], step=Δx[1])
   gridy = range(gridmin[2], gridmax[2], step=Δx[2])
   gridz = range(gridmin[3], gridmax[3], step=Δx[3])

   gridx, gridy, gridz
end

function getinterp(A, gridx, gridy, gridz)
   @assert size(A,1) == 3 && ndims(A) == 4 "Only support 3D force field!"

   Ax = @view A[1,:,:,:]
   Ay = @view A[2,:,:,:]
   Az = @view A[3,:,:,:]

   itp = extrapolate(interpolate(Ax, BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpx = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Ay, BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpy = scale(itp, gridx, gridy, gridz)

   itp = extrapolate(interpolate(Az, BSpline(Cubic(Interpolations.Line(OnGrid())))), NaN)
   interpz = scale(itp, gridx, gridy, gridz)

   # Return field value at a given location.
   function get_field(xu)
      r = @view xu[1:3]

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
         throw(ArgumentError("All methods for the field function $name had too many arguments."))
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
   Field{itd, F}(field_function::F) where {itd, F} = isa(itd, Bool) ? new(field_function) : throw(ArgumentError("itd must be a boolean."))
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

"""
    prepare(grid, E, B; species=Proton) -> (q, m, E, B)

Return a tuple consists of particle charge, mass for a prescribed `species` and interpolated
EM field functions.
"""
function prepare(grid::CartesianGrid, E::TE, B::TB; species::Species=Proton, q=1.0, m=1.0) where {TE, TB}
   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)

   E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz) : Field(E)
   B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz) : Field(B)

   q, m, E, B
end

"""
    prepare(grid, E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, interpolated EM field functions, and external force `F`.
"""
function prepare(grid::CartesianGrid, E::TE, B::TB, F::TF; species::Species=Proton, q=1.0, m=1.0) where {TE, TB, TF}
   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)

   E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz) : Field(E)
   B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz) : Field(B)
   F = TF <: AbstractArray ? getinterp(F, gridx, gridy, gridz) : Field(F)

   q, m, E, B, F
end

"""
    prepare(E, B; species=Proton, q=1.0, m=1.0) -> (q, m, E, B)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m` and analytic EM field functions. Prescribed `species` are `Electron` and
`Proton`; other species can be manually specified with `species=Ion/User`, `q` and `m`.
"""
function prepare(E, B; species::Species=Proton, q=1.0, m=1.0)
   q, m = getchargemass(species, q, m)
   B = Field(B)
   E = Field(E)

   q, m, E, B
end

"""
    prepare(E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, analytic EM field functions, and external force `F`.
"""
function prepare(E, B, F; species::Species=Proton, q=1.0, m=1.0)
   q, m = getchargemass(species, q, m)
   B = Field(B)
   E = Field(E)
   F = Field(F)

   q, m, E, B, F
end

"""
    trace!(dy, y, p::TPTuple, t)

ODE equations for charged particle moving in static EM field with in-place form.
"""
function trace!(dy, y, p::TPTuple, t)
   q, m, E, B = p
   v = @view y[4:6]
   dy[1:3] = v
   dy[4:6] = q/m*(E(y, t) + v × (B(y, t)))
end

"""
    trace(y, p::TPTuple, t) -> SVector{6, Float64}

ODE equations for charged particle moving in static EM field with out-of-place form.
"""
function trace(y, p::TPTuple, t)
   q, m, E, B = p
   v = @view y[4:6]
   dx, dy, dz = v
   dux, duy, duz = q/m*(E(y, t) + v × (B(y, t)))
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

"""
    trace!(dy, y, p::FullTPTuple, t)

ODE equations for charged particle moving in static EM field and external force
field with in-place form.
"""
function trace!(dy, y, p::FullTPTuple, t)
   q, m, E, B, F = p
   v = @view y[4:6]
   dy[1:3] = v
   dy[4:6] = (q*(E(y, t) + v × (B(y, t))) + F(y, t)) / m
end

"""
    trace(y, p::FullTPTuple, t) -> SVector{6, Float64}

ODE equations for charged particle moving in static EM field and external force
field with out-of-place form.
"""
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

ODE equations for relativistic charged particle moving in static EM field with in-place form.
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

ODE equations for relativistic charged particle moving in static EM field with out-of-place form.
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

end