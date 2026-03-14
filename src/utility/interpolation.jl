# Field interpolations.

"""
    AbstractFieldInterpolator

Abstract type for all field interpolators.
"""
abstract type AbstractFieldInterpolator <: Function end

"""
    FieldInterpolator{T}

A callable struct that wraps a 3D interpolation object.
"""
struct FieldInterpolator{T} <: AbstractFieldInterpolator
    itp::T
end

const FieldInterpolator3D = FieldInterpolator

@inline @inbounds function (fi::FieldInterpolator)(xu)
    return fi.itp((xu[1], xu[2], xu[3]))
end

@inline function (fi::FieldInterpolator)(xu, t)
    return fi(xu)
end

@inline @inbounds function jacobian(fi::FieldInterpolator, xu)
    # gradient returns a tuple of partial derivatives:
    # (df/dx, df/dy, df/dz)
    grads = gradient(fi.itp, (xu[1], xu[2], xu[3]))
    # For a vector field f = [f1, f2, f3], each grad is an SVector{3}
    # We want J_ij = ∂fi/∂xj, which means the columns are the partial derivatives.
    return hcat(grads...)
end

Adapt.adapt_structure(to, fi::FieldInterpolator) = FieldInterpolator(Adapt.adapt(to, fi.itp))

"""
    FieldInterpolator2D{T}

A callable struct that wraps a 2D interpolation object.
"""
struct FieldInterpolator2D{T} <: AbstractFieldInterpolator
    itp::T
end

@inline @inbounds function (fi::FieldInterpolator2D)(xu)
    # 2D interpolation usually involves x and y
    return fi.itp((xu[1], xu[2]))
end

@inline function (fi::FieldInterpolator2D)(xu, t)
    return fi(xu)
end

@inline @inbounds function jacobian(fi::FieldInterpolator2D, xu)
    grads = gradient(fi.itp, (xu[1], xu[2]))
    return hcat(grads...)
end

Adapt.adapt_structure(to, fi::FieldInterpolator2D) = FieldInterpolator2D(Adapt.adapt(to, fi.itp))

"""
    FieldInterpolator1D{T}

A callable struct that wraps a 1D interpolation object.
"""
struct FieldInterpolator1D{T} <: AbstractFieldInterpolator
    itp::T
    dir::Int
end

@inline @inbounds function (fi::FieldInterpolator1D)(xu)
    return fi.itp(xu[fi.dir])
end

@inline function (fi::FieldInterpolator1D)(xu, t)
    return fi(xu)
end

@inline @inbounds function jacobian(fi::FieldInterpolator1D, xu)
    # deriv1 returns a viewer that evaluates to the derivative
    return hcat(deriv1(fi.itp)(xu[fi.dir]))
end

Adapt.adapt_structure(to, fi::FieldInterpolator1D) = FieldInterpolator1D(Adapt.adapt(to, fi.itp), fi.dir)

function jacobian(f::Function, x)
    return ForwardDiff.jacobian(f, x)
end

"""
    SphericalFieldInterpolator{T}

A callable struct for spherical grid interpolation (scalar or combined vector).
"""
struct SphericalFieldInterpolator{T} <: AbstractFieldInterpolator
    itp::T
end

@inline function (fi::SphericalFieldInterpolator)(xu)
    rθϕ = cart2sph(xu)
    res = fi.itp(rθϕ)
    if length(res) > 1
        # Convert vector result from spherical to cartesian basis
        Br, Bθ, Bϕ = res
        return @inbounds sph_to_cart_vector(Br, Bθ, Bϕ, rθϕ[2], rθϕ[3])
    else
        return res
    end
end

@inline function (fi::SphericalFieldInterpolator)(xu, t)
    return fi(xu)
end

Adapt.adapt_structure(to, fi::SphericalFieldInterpolator) = SphericalFieldInterpolator(Adapt.adapt(to, fi.itp))

function _get_extrap_mode(bc, T::Type)
    if bc == 2
        return Extrap(:wrap)
    elseif bc == 3
        return Extrap(:clamp)
    else
        if T <: SVector
            return Extrap(:fill; fill_value = SVector{3, eltype(T)}(NaN, NaN, NaN))
        else
            return Extrap(:fill; fill_value = T(NaN))
        end
    end
end

function _fastinterp(grids, A, order, bc)
    extrap_mode = _get_extrap_mode(bc, eltype(A))
    if order == 1
        return linear_interp(grids, A; extrap = extrap_mode)
    elseif order == 2
        return quadratic_interp(grids, A; extrap = extrap_mode)
    elseif order == 3
        return cubic_interp(grids, A; extrap = extrap_mode)
    else
        return constant_interp(grids, A; extrap = extrap_mode)
    end
end


function _fastinterp_spherical(grids, A, order, ϕ_inclusive::Bool)
    # r and θ always extrapolate with NaN, ϕ is always periodic.
    T = eltype(A)
    fill_value = T <: SVector ? SVector{3, eltype(T)}(NaN, NaN, NaN) : T(NaN)
    extrap = (Extrap(:fill; fill_value), Extrap(:fill; fill_value), Extrap(:wrap))
    if order == 1
        return linear_interp(grids, A; extrap)
    elseif order == 2
        # FastInterpolations quadratic does not support PeriodicBC.
        # We use ZeroCurvBC for all axes; Extrap(:wrap) handles periodicity.
        bc_quad = (ZeroCurvBC(), ZeroCurvBC(), ZeroCurvBC())
        return quadratic_interp(grids, A; bc = bc_quad, extrap)
    elseif order == 3
        # When the ϕ grid spans [0, 2π] inclusively, PeriodicBC() (inclusive endpoint) is
        # used for the cubic case; otherwise PeriodicBC(endpoint=:exclusive) is used.
        ϕ_bc = ϕ_inclusive ? PeriodicBC() : PeriodicBC(; endpoint = :exclusive, period = 2π)
        bc_cubic = (ZeroCurvBC(), ZeroCurvBC(), ϕ_bc)
        return cubic_interp(grids, A; bc = bc_cubic, extrap)
    else
        return constant_interp(grids, A; extrap)
    end
end

@inline build_interpolator(A, grid1, args...) = build_interpolator(CartesianGrid, A, grid1, args...)

"""
    build_interpolator(gridtype, A, grids..., order::Int=1, bc::Int=1)
    build_interpolator(A, grids..., order::Int=1, bc::Int=1)

Return a function for interpolating field array `A` on the given grids.

# Arguments

  - `gridtype`: `CartesianGrid`, `RectilinearGrid` or `StructuredGrid`. Usually determined by the number of grids.
  - `A`: field array. For vector field, the first dimension should be 3 if it's not an SVector wrapper.
  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Clamp (flat extrapolation).

# Notes
The input array `A` may be modified in-place for memory optimization.
"""
function build_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(CartesianGrid, As, gridx, gridy, gridz, order, bc)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    itp = _fastinterp((gridx, gridy, gridz), A, order, bc)
    return FieldInterpolator(itp)
end

function build_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(RectilinearGrid, As, gridx, gridy, gridz, order, bc)
end

function build_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 1
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    if order != 1
        throw(ArgumentError("RectilinearGrid (CartesianNonUniform) only supports order=1 (Linear) interpolation."))
    end

    itp = _fastinterp((gridx, gridy, gridz), A, order, bc)
    return FieldInterpolator(itp)
end

function build_interpolator(
        ::Type{<:StructuredGrid}, A::AbstractArray{T, 4},
        gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(StructuredGrid, As, gridr, gridθ, gridϕ, order, bc)
end

function build_interpolator(
        ::Type{<:StructuredGrid}, A::AbstractArray{T, 3},
        gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    r_min, r_max = extrema(gridr)
    θ_min, θ_max = extrema(gridθ)
    ϕ_min, ϕ_max = extrema(gridϕ)

    @assert r_min >= 0 "r must be non-negative."
    @assert θ_min >= 0 && θ_max <= π "θ must be within [0, π]."
    @assert ϕ_min >= 0 && ϕ_max <= 2π "ϕ must be within [0, 2π]."

    has_0 = isapprox(ϕ_min, 0, atol = 1.0e-5)
    has_2pi = isapprox(ϕ_max, 2π, atol = 1.0e-5)
    ϕ_inclusive = has_0 && has_2pi  # grid covers full period with both endpoints

    itp = _fastinterp_spherical((gridr, gridθ, gridϕ), A, order, ϕ_inclusive)
    return SphericalFieldInterpolator(itp)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A,
        gridx::AbstractVector, gridy::AbstractVector, order::Int = 1, bc::Int = 1
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 2 "Inconsistent 2D force field and grid! Expected 2D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 3 "Inconsistent 2D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    itp = _fastinterp((gridx, gridy), As, order, bc)
    return FieldInterpolator2D(itp)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A, gridx::AbstractVector,
        order::Int = 1, bc::Int = 1; dir = 1
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 1 "Inconsistent 1D force field and grid! Expected 1D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 2 "Inconsistent 1D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    itp = _fastinterp(gridx, As, order, bc)

    return FieldInterpolator1D(itp, dir)
end


# Time-dependent field interpolation.

"""
    LazyTimeInterpolator{T, F, L}

A callable struct for handling time-dependent fields with lazy loading and linear time interpolation.

# Fields

- `times::Vector{T}`: Sorted vector of time points.
- `loader::L`: Function `i -> field` that loads the field at index `i`.
- `buffer::Dict{Int, F}`: Cache for loaded fields.
- `lock::ReentrantLock`: Lock for thread safety.
"""
struct LazyTimeInterpolator{T, F, L} <: Function
    times::Vector{T}
    loader::L
    buffer::Dict{Int, F}
    lock::ReentrantLock
end

function LazyTimeInterpolator(times::AbstractVector, loader::Function)
    # Determine the field type by loading the first field
    f1 = loader(1)
    return _LazyTimeInterpolator(times, loader, f1)
end

function _LazyTimeInterpolator(times::AbstractVector, loader::Function, f1::F) where {F}
    buffer = Dict{Int, F}(1 => f1)
    lock = ReentrantLock()
    return LazyTimeInterpolator{eltype(times), F, typeof(loader)}(
        times, loader, buffer, lock
    )
end

function (itp::LazyTimeInterpolator)(x, t)
    # Find the time interval [t1, t2] such that t1 <= t <= t2 (assume times is sorted)
    idx = searchsortedlast(itp.times, t)

    # Handle out-of-bounds
    if idx == 0
        return _get_field!(itp, 1)(x) # clamp to start
    elseif idx >= length(itp.times)
        return _get_field!(itp, length(itp.times))(x) # clamp to end
    end

    t1 = itp.times[idx]
    t2 = itp.times[idx + 1]

    w = (t - t1) / (t2 - t1) # linear weights

    # Load fields (lazily)
    f1 = _get_field!(itp, idx)
    f2 = _get_field!(itp, idx + 1)

    return (1 - w) * f1(x) + w * f2(x)
end

function _get_field!(itp::LazyTimeInterpolator, idx::Int)
    return lock(itp.lock) do
        if !haskey(itp.buffer, idx)
            # Remove far-away indices
            filter!(p -> abs(p.first - idx) <= 1, itp.buffer)

            field = itp.loader(idx)
            itp.buffer[idx] = field
        end
        return itp.buffer[idx]
    end
end

function jacobian(itp::LazyTimeInterpolator, x, t)
    idx = searchsortedlast(itp.times, t)

    # Handle out-of-bounds (clamp)
    if idx <= 0
        return jacobian(_get_field!(itp, 1), x)
    elseif idx >= length(itp.times)
        return jacobian(_get_field!(itp, length(itp.times)), x)
    end

    t1 = itp.times[idx]
    t2 = itp.times[idx + 1]
    w = (t - t1) / (t2 - t1)

    # Loading the fields here so we can get their jacobians
    f1 = _get_field!(itp, idx)
    f2 = _get_field!(itp, idx + 1)

    # Jacobian is linear combination of step jacobians
    return (1 - w) * jacobian(f1, x) + w * jacobian(f2, x)
end

function derivative_t(itp::LazyTimeInterpolator, x, t)
    idx = searchsortedlast(itp.times, t)

    if idx <= 0 || idx >= length(itp.times)
        return zero(itp.buffer[1](x))
    end

    t1 = itp.times[idx]
    t2 = itp.times[idx + 1]

    f1 = _get_field!(itp, idx)
    f2 = _get_field!(itp, idx + 1)

    return (f2(x) - f1(x)) / (t2 - t1)
end
