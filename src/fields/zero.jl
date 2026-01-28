struct ZeroField <: AbstractField{false} end

struct ZeroVector end

# ZeroVector operations
(+)(::ZeroVector, x) = x
(+)(x, ::ZeroVector) = x
(+)(::ZeroVector, ::ZeroVector) = ZeroVector()
(-)(::ZeroVector, x) = -x
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
Base.setindex!(A::AbstractArray, ::ZeroVector, I...) = fill!(view(A, I...), zero(eltype(A)))
Base.setindex!(A::Array{Any}, z::ZeroVector, i::Int) =
    invoke(Base.setindex!, Tuple{Array{Any}, Any, Int}, A, z, i)

function Base.setindex!(A::Array, ::ZeroVector, i::Int)
    return A[i] = zero(eltype(A))
end

Base.getindex(::ZeroVector, I...) = 0
Base.length(::ZeroVector) = 3
Base.iterate(::ZeroVector, state = 1) = state > 3 ? nothing : (0, state + 1)

# Field interface
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
