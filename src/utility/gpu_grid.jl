# Device-compatible grid field interpolators.

"""
    GPUUniformAxis{T, B}

Uniformly spaced grid axis of a [`GPUGrid`](@ref TestParticle.GPUGrid); a cell is located
with a multiply.
"""
struct GPUUniformAxis{T, B}
    x0::T
    inv_dx::T
    n::Int32
    xmax::T
    bc::B
end

function GPUUniformAxis(x0::T, dx::T, n::Integer, bc) where {T}
    return GPUUniformAxis(x0, inv(dx), Int32(n), x0 + T(n - 1) * dx, bc)
end

"""
    GPUNonUniformAxis{T, V<:AbstractVector{T}, B}

Non-uniformly spaced grid axis of a [`GPUGrid`](@ref TestParticle.GPUGrid); a cell is
bracketed with a binary search.
"""
struct GPUNonUniformAxis{T, V <: AbstractVector{T}, B}
    coords::V
    x0::T
    xmax::T
    n::Int32
    bc::B
end

GPUNonUniformAxis(coords::AbstractVector, bc) =
    GPUNonUniformAxis(coords, first(coords), last(coords), Int32(length(coords)), bc)

Adapt.@adapt_structure GPUNonUniformAxis

# `D` holds the spanned location components as a type parameter so that indexing a
# location with them compiles to constant offsets.
struct CartesianCoords{D} end

struct SphericalCoords end

"""
    GPUGrid{N, C, X, V, A<:AbstractArray{V, N}} <: AbstractFieldInterpolator

Multilinear interpolator over an `N`-dimensional grid of `V` values that can be evaluated
inside a device kernel. Each axis is a [`GPUUniformAxis`](@ref TestParticle.GPUUniformAxis)
or a [`GPUNonUniformAxis`](@ref TestParticle.GPUNonUniformAxis) carrying its boundary
condition; queries outside a filling axis return `fill`. With `CartesianCoords{D}`, the grid
spans the location components `D`; with `SphericalCoords`, it spans `(r, θ, ϕ)` and
rotates a vector result back into the Cartesian basis.
"""
struct GPUGrid{N, C, X <: NTuple{N, Any}, V, A <: AbstractArray{V, N}} <: AbstractFieldInterpolator
    data::A
    axes::X
    fill::V
    coords::C
end

Adapt.@adapt_structure GPUGrid

@inline function _nan_fill(::Type{V}) where {V}
    if V <: SVector
        T = eltype(V)
        return SVector{length(V), T}(ntuple(_ -> T(NaN), length(V)))
    elseif V <: Number
        return V(NaN)
    else
        return zero(V)
    end
end

@inline _lerp(c0, c1, w) = muladd(w, c1 - c0, c0)

@inline function _search_bracket(coords::AbstractVector{T}, val::T, n::Int32) where {T}
    c_first = @inbounds coords[1]
    c_last = @inbounds coords[n]
    val <= c_first && return (Int32(1), zero(T))
    val >= c_last && return (n - Int32(1), one(T))

    low = Int32(1)
    high = n
    @inbounds while low <= high
        mid = (low + high) >> 1
        if coords[mid] <= val
            low = mid + Int32(1)
        else
            high = mid - Int32(1)
        end
    end
    i = clamp(high, Int32(1), n - Int32(1))
    @inbounds begin
        c_i = coords[i]
        c_next = coords[i + 1]
    end
    w = (val - c_i) / (c_next - c_i)
    return (i, clamp(w, zero(T), one(T)))
end

@inline function _axis_eval(ax::GPUUniformAxis{T}, val::Real) where {T}
    v = T(val)
    isnan(v) && return (Int32(1), zero(T), true)
    v_adj, out = _apply_extrap(v, ax.x0, ax.xmax, ax.bc)
    out && return (Int32(1), zero(T), true)
    fv = clamp((v_adj - ax.x0) * ax.inv_dx, zero(T), T(ax.n - Int32(1)))
    i0 = min(unsafe_trunc(Int32, fv), ax.n - Int32(2))
    w = fv - T(i0)
    return (i0 + Int32(1), w, false)
end

@inline function _axis_eval(ax::GPUNonUniformAxis{T}, val::Real) where {T}
    v = T(val)
    isnan(v) && return (Int32(1), zero(T), true)
    v_adj, out = _apply_extrap(v, ax.x0, ax.xmax, ax.bc)
    out && return (Int32(1), zero(T), true)
    i, w = _search_bracket(ax.coords, v_adj, ax.n)
    return (i, w, false)
end

# Lerps along the last remaining axis, prepending its corner index to `I`, so that the
# innermost lerps pair values adjacent in memory.
@inline _multilinear(data, ::Tuple{}, ::Tuple{}, I...) = @inbounds data[I...]
@inline function _multilinear(data, i::Tuple, w::Tuple, I...)
    lo = _multilinear(data, Base.front(i), Base.front(w), i[end], I...)
    hi = _multilinear(data, Base.front(i), Base.front(w), i[end] + Int32(1), I...)
    return _lerp(lo, hi, w[end])
end

@inline function _interpolate(g::GPUGrid{N}, q::NTuple{N, Real}) where {N}
    cells = map(_axis_eval, g.axes, q)
    any(c -> c[3], cells) && return g.fill
    return _multilinear(g.data, map(c -> c[1], cells), map(c -> c[2], cells))
end

@inline (g::GPUGrid{N, <:CartesianCoords})(q::NTuple{N, Real}) where {N} = _interpolate(g, q)
@inline (g::GPUGrid{N, CartesianCoords{D}})(xu::AbstractVector) where {N, D} =
    _interpolate(g, map(d -> xu[d], D))

@inline function (g::GPUGrid{3, SphericalCoords})(q::NTuple{3, Real})
    r, θ, ϕ = cart2sph(q...)
    return _to_cartesian(_interpolate(g, (r, θ, ϕ)), θ, ϕ)
end
@inline (g::GPUGrid{3, SphericalCoords})(xu::AbstractVector) = g((xu[1], xu[2], xu[3]))

@inline (g::GPUGrid)(x::Real, xs::Real...) = g((x, xs...))
@inline (g::GPUGrid)(xu::Union{AbstractVector, Tuple}, t) = g(xu)

@inline _to_cartesian(res::SVector{3}, θ, ϕ) = sph2cartvec(res[1], res[2], res[3], θ, ϕ)
@inline _to_cartesian(res, θ, ϕ) = res

_gpu_coords(::FieldInterpolator) = CartesianCoords{(1, 2, 3)}()
_gpu_coords(::FieldInterpolator2D) = CartesianCoords{(1, 2)}()
_gpu_coords(fi::FieldInterpolator1D) = CartesianCoords{(fi.dir,)}()
_gpu_coords(::SphericalFieldInterpolator) = SphericalCoords()

# Equally spaced vectors get a uniform axis: locating a cell is then a multiply, not a search.
function _to_gpu_axis(g::AbstractVector, backend, bc)
    T = float(eltype(g))
    n = length(g)
    x0 = T(first(g))
    dx = g isa AbstractRange ? T(step(g)) : (T(last(g)) - x0) / (n - 1)
    uniform = g isa AbstractRange ||
        all(i -> abs(g[i] - (x0 + (i - 1) * dx)) <= sqrt(eps(T)) * abs(dx), 1:n)
    uniform && return GPUUniformAxis(x0, dx, n, bc)
    coords = collect(T, g)
    return GPUNonUniformAxis(Adapt.adapt(backend, coords), first(coords), last(coords), Int32(n), bc)
end

function _to_gpu_grid(
        fi::Union{
            FieldInterpolator, FieldInterpolator2D, FieldInterpolator1D,
            SphericalFieldInterpolator,
        },
        backend
    )
    bcs = _device_bcs(fi.spec)
    axes = map((g, bc) -> _to_gpu_axis(g, backend, bc), fi.spec.grids, bcs)
    data = Adapt.adapt(backend, _itp_data(fi.itp))
    return GPUGrid(data, axes, _fill_value(bcs, eltype(data)), _gpu_coords(fi))
end
