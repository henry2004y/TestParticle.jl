# Construction of tracing parameters.

"""
Judge whether the field function is time dependent.
"""
is_time_dependent(f::Function) = applicable(f, zeros(6), 0.0) ? true : false

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
   function Field{itd, F}(field_function::F) where {itd, F}
      isa(itd, Bool) ? new(field_function) : throw(ArgumentError("itd must be a boolean."))
   end
end

Field(f::Function) = Field{is_time_dependent(f), typeof(f)}(f)

(f::AbstractField{true})(xu, t) = f.field_function(xu, t)
function (f::AbstractField{true})(xu)
   throw(ArgumentError("Time-dependent field function must have a time argument."))
end
(f::AbstractField{false})(xu, t) = SVector{3}(f.field_function(xu))
(f::AbstractField{false})(xu) = SVector{3}(f.field_function(xu))

function Base.show(io::IO, f::Field)
   println(io, "Field with interpolation support")
   println(io, "Time-dependent: ", is_time_dependent(f.field_function))
end


"""
The type of parameter tuple for guiding center problem.
"""
GCTuple = Tuple{Float64, Float64, Float64, AbstractField, AbstractField}

get_q2m(param) = param[1]
get_BField(param) = param[4]
get_EField(param) = param[3]

prepare_field(f, x...; kwargs...) = Field(f)
function prepare_field(f::AbstractArray, x...; order, bc, kw...)
   Field(getinterp(f, x..., order, bc; kw...))
end

function _prepare(E, B, F, args...; species::Species = Proton, q = 1.0, m = 1.0, kw...)
   q, m = getchargemass(species, q, m)
   q2m = q / m
   fE = prepare_field(E, args...; kw...)
   fB = prepare_field(B, args...; kw...)
   fF = prepare_field(F, args...; kw...)
   q2m, m, fE, fB, fF
end

"""
   prepare(args...; kwargs...) -> (q2m, m, E, B, F)
   prepare(E, B, F = ZeroField(); kwargs...)
   prepare(grid::CartesianGrid, E, B, F = ZeroField(); kwargs...)
   prepare(x, E, B, F = ZeroField(); dir = 1, kwargs...)
   prepare(x, y, E, B, F = ZeroField(); kwargs...)
   prepare(x, y, z, E, B, F = ZeroField(); kwargs...)

Return a tuple consists of particle charge-mass ratio for a prescribed `species` of charge `q` and mass `m`,
mass `m` for a prescribed `species`, analytic/interpolated EM field functions, and external force `F`.

Prescribed `species` are `Electron` and `Proton`;
other species can be manually specified with `species=Ion/User`, `q` and `m`.

Direct range input for uniform grid in 1/2/3D is supported.
For 1D grid, an additional keyword `dir` is used for specifying the spatial direction, 1 -> x, 2 -> y, 3 -> z.

# Keywords

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic.
  - `species::Species=Proton`: particle species.
  - `q=1.0`: particle charge. Only works when `Species=User`.
  - `m=1.0`: particle mass. Only works when `Species=User`.
"""
function prepare(grid::CartesianGrid, E, B, F = ZeroField(); order = 1, bc = 1, kw...)
   _prepare(E, B, F, makegrid(grid)...; order, bc, kw...)
end

function prepare(x::T, y::T, E, B, F = ZeroField(); order = 1,
      bc = 1, kw...) where {T <: AbstractRange}
   _prepare(E, B, F, x, y; order, bc, kw...)
end

function prepare(x::T, y::T, z::T, E, B, F = ZeroField(); order = 1,
      bc = 1, kw...) where {T <: AbstractRange}
   _prepare(E, B, F, x, y, z; order, bc, kw...)
end

function prepare(x::AbstractRange, E, B, F = ZeroField(); order = 1, bc = 3, dir = 1, kw...)
   _prepare(E, B, F, x; order, bc, dir, kw...)
end
prepare(E, B, F = ZeroField(); kw...) = _prepare(E, B, F; kw...)

struct ZeroField <: AbstractField{false} end

struct ZeroVector end

# ZeroVector operations
(+)(::ZeroVector, x) = x
(+)(x, ::ZeroVector) = x
(+)(::ZeroVector, ::ZeroVector) = ZeroVector()
(*)(::ZeroVector, _) = ZeroVector()
(*)(_, ::ZeroVector) = ZeroVector()
(/)(::ZeroVector, _) = ZeroVector()
(×)(::ZeroVector, _) = ZeroVector()
(×)(_, ::ZeroVector) = ZeroVector()

# Convert ZeroVector to SVector{3} when needed
(::Type{T})(::ZeroVector) where {T <: StaticArray} = T(0, 0, 0)

# Make ZeroVector work with array assignment
Base.setindex!(A::AbstractArray, ::ZeroVector, I...) = fill!(view(A, I...), 0)
Base.getindex(::ZeroVector, I...) = 0

# Field interface
Field(x::ZeroField) = x
(::ZeroField)(y, t) = ZeroVector()
(::ZeroField)(_) = ZeroVector()
