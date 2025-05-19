# Construction of tracing parameters.

"Judge whether the field function is time dependent."
function is_time_dependent(f::Function)
	applicable(f, zeros(6), 0.0) ? true : false
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
(f::AbstractField{true})(xu) =
	throw(ArgumentError("Time-dependent field function must have a time argument."))
(f::AbstractField{false})(xu, t) = SVector{3}(f.field_function(xu))
(f::AbstractField{false})(xu) = SVector{3}(f.field_function(xu))

function Base.show(io::IO, f::Field)
	println(io, "Field with interpolation support")
	println(io, "Time-dependent: ", is_time_dependent(f.field_function))
end

"The type of parameter tuple for full test particle problem."
FullTPTuple = Tuple{Float64, Float64, AbstractField, AbstractField, AbstractField}

"The type of parameter tuple for normal test particle problem."
TPTuple = Tuple{Float64, AbstractField, AbstractField}

"The type of parameter tuple for normalized test particle problem."
TPNormalizedTuple = Tuple{AbstractFloat, AbstractField, AbstractField}

"The type of parameter tuple for guiding center problem."
GCTuple = Tuple{Float64, Float64, Float64, AbstractField, AbstractField}


"""
	 prepare(grid::CartesianGrid, E, B; kwargs...) -> (q2m, E, B)

Return a tuple consists of particle charge-mass ratio for a prescribed `species` and
interpolated EM field functions.

# Keywords
- `order::Int=1`: order of interpolation in [1,2,3].
- `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic.
- `species::Species=Proton`: particle species.
- `q::AbstractFloat=1.0`: particle charge. Only works when `Species=User`.
- `m::AbstractFloat=1.0`: particle mass. Only works when `Species=User`.

	 prepare(grid::CartesianGrid, E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, interpolated EM field functions, and external force `F`.

	 prepare(x::AbstractRange, y::AbstractRange, z::AbstractRange, E, B; kwargs...) -> (q2m, E, B)
	 prepare(x, y, E, B; kwargs...) -> (q2m, E, B)

	 prepare(x::AbstractRange, E, B; kwargs...) -> (q2m, E, B)

1D grid. An additional keyword `dir` is used for specifying the spatial direction, 1 -> x, 2 -> y, 3 -> z.

Direct range input for uniform grid in 2/3D is also accepted.

	 prepare(E, B; kwargs...) -> (q2m, E, B)

Return a tuple consists of particle charge-mass ratio for a prescribed `species` of charge
`q` and mass `m` and analytic EM field functions. Prescribed `species` are `Electron` and
`Proton`; other species can be manually specified with `species=Ion/User`, `q` and `m`.

	 prepare(E, B, F; kwargs...) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, analytic EM field functions, and external force `F`.
"""
function prepare(grid::CartesianGrid, E::TE, B::TB; species::Species = Proton,
	q::AbstractFloat = 1.0, m::AbstractFloat = 1.0, order::Int = 1, bc::Int = 1,
) where {TE, TB}

	q, m = getchargemass(species, q, m)

	if paramdim(grid) == 3
		gridx, gridy, gridz = makegrid(grid)
		E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz, order, bc) : E
		B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz, order, bc) : B
	elseif paramdim(grid) == 2
		gridx, gridy = makegrid(grid)
		E = TE <: AbstractArray ? getinterp(E, gridx, gridy, order, bc) : E
		B = TB <: AbstractArray ? getinterp(B, gridx, gridy, order, bc) : B
	end

	q/m, Field(E), Field(B)
end

function prepare(grid::CartesianGrid, E::TE, B::TB, F::TF; species::Species = Proton,
	q::AbstractFloat = 1.0, m::AbstractFloat = 1.0, order::Int = 1, bc::Int = 1,
) where {TE, TB, TF}

	q, m = getchargemass(species, q, m)

	gridx, gridy, gridz = makegrid(grid)

	E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz, order, bc) : E
	B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz, order, bc) : B
	F = TF <: AbstractArray ? getinterp(F, gridx, gridy, gridz, order, bc) : F

	q, m, Field(E), Field(B), Field(F)
end

function prepare(x::T, y::T, E::TE, B::TB; species::Species = Proton,
	q::AbstractFloat = 1.0,
	m::AbstractFloat = 1.0, order::Int = 1, bc::Int = 1) where {T <: AbstractRange, TE, TB}

	q, m = getchargemass(species, q, m)

	E = TE <: AbstractArray ? getinterp(E, x, y, order, bc) : E
	B = TB <: AbstractArray ? getinterp(B, x, y, order, bc) : B

	q/m, Field(E), Field(B)
end

function prepare(x::T, E::TE, B::TB; species::Species = Proton, q::AbstractFloat = 1.0,
	m::AbstractFloat = 1.0, order::Int = 1, bc::Int = 3, dir = 1,
) where {T <: AbstractRange, TE, TB}

	q, m = getchargemass(species, q, m)

	E = TE <: AbstractArray ? getinterp(E, x, order, bc; dir) : E
	B = TB <: AbstractArray ? getinterp(B, x, order, bc; dir) : B

	q/m, Field(E), Field(B)
end

function prepare(x::T, y::T, z::T, E::TE, B::TB;
	species::Species = Proton, q::AbstractFloat = 1.0, m::AbstractFloat = 1.0,
	order::Int = 1,
	bc::Int = 1) where {T <: AbstractRange, TE, TB}

	q, m = getchargemass(species, q, m)

	E = TE <: AbstractArray ? getinterp(E, x, y, z, order, bc) : E
	B = TB <: AbstractArray ? getinterp(B, x, y, z, order, bc) : B

	q/m, Field(E), Field(B)
end

function prepare(
	E,
	B;
	species::Species = Proton,
	q::AbstractFloat = 1.0,
	m::AbstractFloat = 1.0,
)
	q, m = getchargemass(species, q, m)

	q/m, Field(E), Field(B)
end

function prepare(E, B, F; species::Species = Proton, q::AbstractFloat = 1.0,
	m::AbstractFloat = 1.0)
	q, m = getchargemass(species, q, m)

	q, m, Field(E), Field(B), Field(F)
end


import Base: (+), (*), (/), setindex!
import LinearAlgebra: ×
import StaticArrays: StaticArray

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

# Field interface
Field(x::ZeroField) = x
(::ZeroField)(y, t) = ZeroVector()