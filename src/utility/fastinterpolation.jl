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

@inbounds function (fi::FieldInterpolator)(xu)
    return fi.itp((xu[1], xu[2], xu[3]))
end

function (fi::FieldInterpolator)(xu, t)
    return fi(xu)
end

Adapt.adapt_structure(to, fi::FieldInterpolator) = FieldInterpolator(Adapt.adapt(to, fi.itp))

"""
    FieldInterpolator2D{T}

A callable struct that wraps a 2D interpolation object.
"""
struct FieldInterpolator2D{T} <: AbstractFieldInterpolator
    itp::T
end

@inbounds function (fi::FieldInterpolator2D)(xu)
    return fi.itp((xu[1], xu[2]))
end

function (fi::FieldInterpolator2D)(xu, t)
    return fi(xu)
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

@inbounds function (fi::FieldInterpolator1D)(xu)
    return fi.itp((xu[fi.dir],))
end

function (fi::FieldInterpolator1D)(xu, t)
    return fi(xu)
end

Adapt.adapt_structure(to, fi::FieldInterpolator1D) = FieldInterpolator1D(Adapt.adapt(to, fi.itp), fi.dir)

"""
    SphericalFieldInterpolator{T}

A callable struct for spherical grid interpolation (scalar or combined vector).
"""
struct SphericalFieldInterpolator{T} <: AbstractFieldInterpolator
    itp::T
end

function (fi::SphericalFieldInterpolator)(xu)
    rθϕ = cart2sph(xu)
    res = fi.itp(rθϕ)
    if length(res) > 1
        Br, Bθ, Bϕ = res
        return sph_to_cart_vector(Br, Bθ, Bϕ, rθϕ[2], rθϕ[3])
    else
        return res
    end
end

function (fi::SphericalFieldInterpolator)(xu, t)
    return fi(xu)
end

Adapt.adapt_structure(to, fi::SphericalFieldInterpolator) = SphericalFieldInterpolator(Adapt.adapt(to, fi.itp))

function _get_extrap_mode(bc)
    if bc == 2
        return Extrap(:wrap)
    elseif bc == 3
        return Extrap(:constant)
    else
        # TODO: bc == 1 (NaN outside domain) requires native FastInterpolations.jl support.
        # Once available, replace this with the appropriate Extrap mode and remove the
        # manual bounds-checking code that was previously in the FieldInterpolator call methods.
        return NoExtrap()
    end
end

function _fastinterp(grids, A, order, bc)
    extrap_mode = _get_extrap_mode(bc)
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

@inline build_interpolator(A, grid1, args...) =
    build_interpolator(FastInterpolationsBackend(), CartesianGrid, A, grid1, args...)
@inline build_interpolator(gridtype::Type, A, args...) =
    build_interpolator(FastInterpolationsBackend(), gridtype, A, args...)

"""
    build_interpolator(backend, gridtype, A, grids..., order::Int=1, bc::Int=3)
    build_interpolator(gridtype, A, grids..., order::Int=1, bc::Int=3)
    build_interpolator(A, grids..., order::Int=1, bc::Int=3)

Return a function for interpolating field array `A` on the given grids.

# Arguments

  - `backend`: `FastInterpolationsBackend()` (default) or `InterpolationsBackend()`.
  - `gridtype`: `CartesianGrid`, `RectilinearGrid` or `StructuredGrid`. Usually determined by the number of grids.
  - `A`: field array. For vector field, the first dimension should be 3 if it's not an SVector wrapper.
  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=3`: type of boundary conditions, 1 -> NaN (not yet native; requires FastInterpolations support), 2 -> periodic, 3 -> Flat.

# Notes
The input array `A` may be modified in-place for memory optimization.
"""
function build_interpolator(
        b::FastInterpolationsBackend, ::Type{<:CartesianGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 3
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(b, CartesianGrid, As, gridx, gridy, gridz, order, bc)
end

function build_interpolator(
        ::FastInterpolationsBackend, ::Type{<:CartesianGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 3
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    itp = _fastinterp((gridx, gridy, gridz), A, order, bc)
    return FieldInterpolator(itp)
end

function build_interpolator(
        b::FastInterpolationsBackend, ::Type{<:RectilinearGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 3
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(b, RectilinearGrid, As, gridx, gridy, gridz, order, bc)
end

function build_interpolator(
        ::FastInterpolationsBackend, ::Type{<:RectilinearGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc::Int = 3
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
        b::FastInterpolationsBackend, ::Type{<:StructuredGrid}, A::AbstractArray{T, 4},
        gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 3
    ) where {T}
    @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(b, StructuredGrid, As, gridr, gridθ, gridϕ, order, bc)
end

function build_interpolator(
        ::FastInterpolationsBackend, ::Type{<:StructuredGrid}, A::AbstractArray{T, 3},
        gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 3
    ) where {T}
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    end
    is_uniform_r = gridr isa AbstractRange && !(gridr isa Base.LogRange)

    if is_uniform_r
        itp = _fastinterp((gridr, gridθ, gridϕ), A, order, bc)
    else
        gridϕ, A = _ensure_full_phi(gridϕ, A)
        itp = _fastinterp((gridr, gridθ, gridϕ), A, order, bc)
    end

    return SphericalFieldInterpolator(itp)
end

function build_interpolator(
        ::FastInterpolationsBackend, ::Type{<:CartesianGrid}, A,
        gridx::AbstractVector, gridy::AbstractVector, order::Int = 1, bc::Int = 3
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
        ::FastInterpolationsBackend, ::Type{<:CartesianGrid}, A, gridx::AbstractVector,
        order::Int = 1, bc::Int = 3; dir = 1
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 1 "Inconsistent 1D force field and grid! Expected 1D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 2 "Inconsistent 1D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    itp = _fastinterp((gridx,), As, order, bc)
    return FieldInterpolator1D(itp, dir)
end

function _ensure_full_phi(gridϕ, A::AbstractArray{T, N}) where {T, N}
    min_phi, max_phi = extrema(gridϕ)
    needs_0 = !isapprox(min_phi, 0, atol = 1.0e-5)
    needs_2pi = !isapprox(max_phi, 2π, atol = 1.0e-5)

    if !needs_0 && !needs_2pi
        return gridϕ, A
    end

    new_grid_vec = collect(gridϕ)
    if needs_0
        pushfirst!(new_grid_vec, 0.0)
    end
    if needs_2pi
        push!(new_grid_vec, 2π)
    end

    phi_dim = N
    new_A = Array{T, N}(undef, (size(A)[1:(end - 1)]..., length(new_grid_vec)))

    start_idx = needs_0 ? 2 : 1
    end_idx = start_idx + length(gridϕ) - 1

    selectdim(new_A, phi_dim, start_idx:end_idx) .= A

    if needs_0
        src_idx_for_0 = size(A, phi_dim)
        selectdim(new_A, phi_dim, 1) .= selectdim(A, phi_dim, src_idx_for_0)
    end

    if needs_2pi
        if needs_0
            selectdim(new_A, phi_dim, size(new_A, phi_dim)) .= selectdim(new_A, phi_dim, 1)
        else # needs 2π only
            selectdim(new_A, phi_dim, size(new_A, phi_dim)) .= selectdim(A, phi_dim, 1)
        end
    end

    return new_grid_vec, new_A
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
