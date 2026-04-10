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
    # partial derivatives: (df/dx, df/dy, df/dz)
    grads = gradient(fi.itp, (xu[1], xu[2], xu[3]))
    # For a vector field f = [f1, f2, f3], each grad is an SVector{3}
    return hcat(grads...) # J_ij = ∂fi/∂xj
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
        # Convert vector result from spherical to Cartesian basis
        Br, Bθ, Bϕ = res
        return @inbounds sph2cartvec(Br, Bθ, Bϕ, rθϕ[2], rθϕ[3])
    else
        return res
    end
end

@inline function (fi::SphericalFieldInterpolator)(xu, t)
    return fi(xu)
end

Adapt.adapt_structure(to, fi::SphericalFieldInterpolator) = SphericalFieldInterpolator(Adapt.adapt(to, fi.itp))

function _fastinterp(grids, A, order, extrap::AbstractExtrap, coeffs = OnTheFly())
    # Ensure FillExtrap value matches eltype(A) for SVector types
    if extrap isa FillExtrap && extrap.fill_value isa Number && isnan(extrap.fill_value)
        T = eltype(A)
        if T <: SVector
            extrap = FillExtrap(SVector{3, eltype(T)}(NaN, NaN, NaN))
        else
            extrap = FillExtrap(T(NaN))
        end
    end

    if order == 1
        return linear_interp(grids, A; extrap)
    elseif order == 3
        return cardinal_interp(grids, A; coeffs, extrap)
    elseif order == 0
        return constant_interp(grids, A; extrap)
    else
        throw(ArgumentError("Interpolation order $order is not supported. Supported orders are 0, 1, and 3."))
    end
end


function _fastinterp_spherical(grids, A, order, coeffs = PreCompute())
    # r and θ always extrapolate with NaN, ϕ is always periodic.
    T = eltype(A)
    fill_value = T <: SVector ? SVector{3, eltype(T)}(NaN, NaN, NaN) : T(NaN)
    extrap = (Extrap(:fill; fill_value), Extrap(:fill; fill_value), Extrap(:wrap))
    if order == 1
        return linear_interp(grids, A; extrap)
    elseif order == 3
        return cardinal_interp(grids, A; coeffs, extrap)
    elseif order == 0
        return constant_interp(grids, A; extrap)
    else
        throw(ArgumentError("Interpolation order $order is not supported. Supported orders are 0, 1, and 3."))
    end
end

function _check_interpolation_consistency(A, grids, order)
    order < 3 && return

    T = eltype(A)
    # Get the underlying element type if it's an SVector
    Tv = T <: SVector ? eltype(T) : T

    # Get the promoted type of coordinates from all grids
    Td = promote_type(map(eltype, grids)...)

    # FastInterpolations ND kernels for non-scalar values (like SVector)
    # strictly require value types to match the promoted grid type.
    if ndims(A) > 1 && T <: SVector && Tv != Td
        throw(
            ArgumentError(
                "High-order interpolation (order >= 3) in FastInterpolations requires " *
                    "field data type to match the promoted type of grid coordinates. " *
                    "Found data type $Tv and grid type $Td. " *
                    "Please convert your field data to $Td or your grids to $Tv."
            )
        )
    end
    return
end

@inline build_interpolator(A::AbstractArray, grid1, args...; kwargs...) = build_interpolator(CartesianGrid, A, grid1, args...; kwargs...)

raw"""
    build_interpolator(gridtype, A, grids..., order::Int=1, bc=FillExtrap(NaN))
    build_interpolator(A, grids..., order::Int=1, bc=FillExtrap(NaN))

Return a function for interpolating field array `A` on the given grids.

# Arguments

  - `gridtype`: `CartesianGrid`, `RectilinearGrid` or `StructuredGrid`. Usually determined by the number of grids.
  - `A`: field array. For vector field, the first dimension should be 3 if it's not an SVector wrapper.
  - `order::Int=1`: order of interpolation in [0,1,3].
  - `bc=FillExtrap(NaN)`: boundary condition type from `FastInterpolations.jl`.
    - `FillExtrap(NaN)`: Fill with NaN (default).
    - `ClampExtrap()`: Clamp (flat extrapolation).
    - `WrapExtrap()`: Exclusive periodic wrapping ($L = N \Delta x$).
  - `coeffs=OnTheFly()`: coefficient strategy for cubic interpolation (order=3). Default is `OnTheFly()`.

# Notes
- The input array `A` may be modified in-place for memory optimization.
"""
function build_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN); coeffs = OnTheFly()
    ) where {T}
    @assert size(A, 1) == 3 "Incompatible 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(CartesianGrid, As, gridx, gridy, gridz, order, bc; coeffs)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN); coeffs = OnTheFly()
    ) where {T}
    _check_interpolation_consistency(A, (gridx, gridy, gridz), order)
    itp = _fastinterp((gridx, gridy, gridz), A, order, bc, coeffs)
    return FieldInterpolator(itp)
end

function build_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN)
    ) where {T}
    @assert size(A, 1) == 3 "Incompatible 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(RectilinearGrid, As, gridx, gridy, gridz, order, bc)
end

function build_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN)
    ) where {T}
    if order != 1
        throw(ArgumentError("RectilinearGrid (CartesianNonUniform) only supports order=1 (Linear) interpolation."))
    end

    itp = _fastinterp((gridx, gridy, gridz), A, order, bc)
    return FieldInterpolator(itp)
end

function build_interpolator(
        ::Type{<:StructuredGrid}, A::AbstractArray{T, 4},
        gridr, gridθ, gridϕ, order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN); coeffs = OnTheFly()
    ) where {T}
    @assert size(A, 1) == 3 "Incompatible 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(StructuredGrid, As, gridr, gridθ, gridϕ, order, bc; coeffs)
end

function build_interpolator(
        ::Type{<:StructuredGrid}, A::AbstractArray{T, 3},
        gridr, gridθ, gridϕ, order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN); coeffs = OnTheFly()
    ) where {T}
    r_min = minimum(gridr)
    θ_min, θ_max = extrema(gridθ)
    ϕ_min, ϕ_max = extrema(gridϕ)

    @assert r_min >= 0 "r must be non-negative."
    @assert θ_min >= 0 && θ_max <= π "θ must be within [0, π]."
    @assert ϕ_min >= 0 && ϕ_max <= 2π "ϕ must be within [0, 2π]."

    _check_interpolation_consistency(A, (gridr, gridθ, gridϕ), order)
    itp = _fastinterp_spherical((gridr, gridθ, gridϕ), A, order, coeffs)
    return SphericalFieldInterpolator(itp)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A,
        gridx::AbstractVector, gridy::AbstractVector, order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN);
        coeffs = OnTheFly()
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 2 "Incompatible 2D force field and grid! Expected 2D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 3 "Incompatible 2D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    _check_interpolation_consistency(As, (gridx, gridy), order)
    itp = _fastinterp((gridx, gridy), As, order, bc, coeffs)
    return FieldInterpolator2D(itp)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A, gridx::AbstractVector,
        order::Int = 1, bc::AbstractExtrap = FillExtrap(NaN); dir = 1, coeffs = OnTheFly()
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 1 "Incompatible 1D force field and grid! Expected 1D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 2 "Incompatible 1D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    _check_interpolation_consistency(As, (gridx,), order)
    itp = _fastinterp(gridx, As, order, bc, coeffs)

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
        return zero(_get_field!(itp, 1)(x))
    end

    t1 = itp.times[idx]
    t2 = itp.times[idx + 1]

    f1 = _get_field!(itp, idx)
    f2 = _get_field!(itp, idx + 1)

    return (f2(x) - f1(x)) / (t2 - t1)
end
