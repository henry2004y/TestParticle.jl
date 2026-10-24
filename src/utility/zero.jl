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

Base.getindex(::ZeroVector, I...) = 0
Base.length(::ZeroVector) = 3
Base.iterate(::ZeroVector, state = 1) = state > 3 ? nothing : (0, state + 1)

# Field interface. Field values are also stored, differentiated and combined with
# `⋅`/`×`/`vcat`, so floats get a real zero vector.
function (::ZeroField)(y)
    T = eltype(y)
    if T <: AbstractFloat || T <: Complex{<:AbstractFloat}
        return zero(SVector{3, T})
    else
        return ZeroVector()
    end
end

(z::ZeroField)(y, t) = z(y)

# Field value as `SVector{3, T}`, for the equations that add up field terms. A `ZeroField`
# term becomes a `ZeroVector`, which the arithmetic drops by dispatch.
@inline _eval_field(::Type{T}, f, x, t) where {T} = SVector{3, T}(f(x, t))
_eval_field(::Type, ::ZeroField, x, t) = ZeroVector()
