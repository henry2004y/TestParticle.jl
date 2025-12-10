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
Base.length(::ZeroVector) = 3
Base.iterate(::ZeroVector, state = 1) = state > 3 ? nothing : (0, state + 1)

# Field interface
Field(x::ZeroField) = x
(::ZeroField)(y, t) = ZeroVector()
(::ZeroField)(_) = ZeroVector()
