struct ZeroField <: AbstractField{false} end

struct ZeroVector end

# ZeroVector operations
(+)(::ZeroVector, x) = x
(+)(x, ::ZeroVector) = x
(+)(::ZeroVector, ::ZeroVector) = ZeroVector()
(-)(::ZeroVector, x) = x
(-)(x, ::ZeroVector) = x
(-)(::ZeroVector, ::ZeroVector) = ZeroVector()
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
function (::ZeroField)(y, t)
    T = eltype(y)
    if T <: AbstractFloat || T <: Complex{<:AbstractFloat}
        return zero(SVector{3, T})
    else
        return ZeroVector()
    end
end

function (::ZeroField)(y)
    T = eltype(y)
    if T <: AbstractFloat || T <: Complex{<:AbstractFloat}
        return zero(SVector{3, T})
    else
        return ZeroVector()
    end
end
