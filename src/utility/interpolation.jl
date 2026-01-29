# Field interpolations.

@inline getinterp(A, grid1, args...) = getinterp(CartesianGrid, A, grid1, args...)

"""
    FieldInterpolator{T}

A callable struct that wraps an interpolation object.
It enables compatibility with `boris_kernel` by exposing the inner `itp` object for GPU adaptation.
"""
struct FieldInterpolator{T} <: Function
    itp::T
end

function (fi::FieldInterpolator)(xu)
    return fi.itp(xu[1], xu[2], xu[3])
end

function getinterp_scalar(A, grid1, grid2, grid3, args...)
    return getinterp_scalar(CartesianGrid, A, grid1, grid2, grid3, args...)
end

"""
    getinterp(::Type{<:CartesianGrid}, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)

Return a function for interpolating field array `A` on the grid given by `gridx`, `gridy`, and `gridz`.

# Arguments

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.
  - `dir::Int`: 1/2/3, representing x/y/z direction.

# Notes
The input array `A` may be modified in-place for memory optimization.
"""
function getinterp(
        ::Type{<:CartesianGrid}, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    else
        @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    end
    return get_interpolator(CartesianGrid, A, gridx, gridy, gridz, order, bc)
end

function getinterp(
        ::Type{<:RectilinearGrid}, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    else
        @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    end
    return get_interpolator(RectilinearGrid, A, gridx, gridy, gridz, order, bc)
end

function getinterp(
        ::Type{<:StructuredGrid}, A, gridr, gridθ, gridϕ,
        order::Int = 1, bc::Int = 3
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 3 "Inconsistent 3D force field and grid! Expected 3D array of SVectors."
    else
        @assert size(A, 1) == 3 && ndims(A) == 4 "Inconsistent 3D force field and grid!"
    end
    return get_interpolator(StructuredGrid, A, gridr, gridθ, gridϕ, order, bc)
end

function getinterp(::Type{<:CartesianGrid}, A, gridx, gridy, order::Int = 1, bc::Int = 2)
    if eltype(A) <: SVector
        @assert ndims(A) == 2 "Inconsistent 2D force field and grid! Expected 2D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 3 "Inconsistent 2D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    itp = _get_interp_object(As, order, bc)

    interp = scale(itp, gridx, gridy)

    # Return field value at a given location.
    function get_field(xu)
        r = @view xu[1:2]

        return interp(r...)
    end

    return get_field
end

function get_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 4},
        gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    ) where {T}
    As = reinterpret(reshape, SVector{3, T}, A)
    return get_interpolator(RectilinearGrid, As, gridx, gridy, gridz, order, bc)
end

function get_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 3},
        gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    ) where {T}
    if order != 1
        throw(ArgumentError("RectilinearGrid (CartesianNonUniform) only supports order=1 (Linear) interpolation."))
    end

    bctype = if bc == 1
        if T <: SVector
            fill(eltype(T)(NaN), 3)
        else
            T(NaN)
        end
    elseif bc == 2
        Periodic()
    else
        Flat()
    end

    itp = extrapolate(interpolate!((gridx, gridy, gridz), A, Gridded(Linear())), bctype)

    function get_field(xu)
        r = @view xu[1:3]
        return itp(r...)
    end

    return get_field
end

function getinterp(::Type{<:CartesianGrid}, A, gridx, order::Int = 1, bc::Int = 3; dir = 1)
    if eltype(A) <: SVector
        @assert ndims(A) == 1 "Inconsistent 1D force field and grid! Expected 1D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 2 "Inconsistent 1D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    itp = _get_interp_object(As, order, bc)

    interp = scale(itp, gridx)

    # Return field value at a given location.
    function get_field(xu)
        r = xu[dir]

        return interp(r)
    end

    return get_field
end

function _get_bspline(order::Int, periodic::Bool)
    gt = OnCell()

    interp_type = if order == 1
        Linear
    elseif order == 2
        Quadratic
    elseif order == 3
        Cubic
    else
        throw(ArgumentError("Unsupported interpolation order!"))
    end

    if periodic
        return BSpline(interp_type(Periodic(gt)))
    else
        # Linear() is special as it doesn't take an argument.
        if interp_type == Linear
            return BSpline(Linear())
        else
            return BSpline(interp_type(Flat(gt)))
        end
    end
end

function _getinterp(Ax, Ay, Az, order::Int, bc::Int)
    itpx = _get_interp_object(Ax, order, bc)
    itpy = _get_interp_object(Ay, order, bc)
    itpz = _get_interp_object(Az, order, bc)

    return itpx, itpy, itpz
end

function _getinterp(gridtype::Type{<:StructuredGrid}, Ax, Ay, Az, order::Int, bc::Int)
    itpr = _get_interp_object(gridtype, Ax, order, bc)
    itpθ = _get_interp_object(gridtype, Ay, order, bc)
    itpϕ = _get_interp_object(gridtype, Az, order, bc)

    return itpr, itpθ, itpϕ
end

"""
    getinterp_scalar(::Type{<:CartesianGrid}, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)

Return a function for interpolating scalar array `A` on the grid given by `gridx`, `gridy`, and `gridz`. Currently only 3D arrays are supported.

# Arguments

  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.
  - `dir::Int`: 1/2/3, representing x/y/z direction.
"""
function getinterp_scalar(
        ::Type{<:CartesianGrid}, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    )
    return get_interpolator(CartesianGrid, A, gridx, gridy, gridz, order, bc)
end

function getinterp_scalar(
        ::Type{<:RectilinearGrid}, A, gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    )
    return get_interpolator(RectilinearGrid, A, gridx, gridy, gridz, order, bc)
end

function getinterp_scalar(
        ::Type{<:StructuredGrid}, A,
        gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 3
    )
    return get_interpolator(StructuredGrid, A, gridr, gridθ, gridϕ, order, bc)
end

"""
    get_interpolator(A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
    get_interpolator(gridtype, A, grid1, grid2, grid3, order::Int=1, bc::Int=1)

Return a function for interpolating field array `A` on the grid given by `gridx`, `gridy`, and `gridz`.

# Arguments

  - `gridtype`: `CartesianGrid`, `RectilinearGrid` or `StructuredGrid`.
  - `A`: field array. For vector field, the first dimension should be 3.
  - `order::Int=1`: order of interpolation in [1,2,3].
  - `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic, 3 -> Flat.

# Notes
The input array `A` may be modified in-place for memory optimization.
"""
function get_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 4},
        gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    ) where {T}
    As = reinterpret(reshape, SVector{3, T}, A)
    return get_interpolator(CartesianGrid, As, gridx, gridy, gridz, order, bc)
end

function get_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 3},
        gridx, gridy, gridz, order::Int = 1, bc::Int = 1
    ) where {T}
    itp = _get_interp_object(A, order, bc)

    interp = scale(itp, gridx, gridy, gridz)

    # Return field value at a given location.
    return FieldInterpolator(interp)
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
            selectdim(new_A, phi_dim, size(new_A, phi_dim)) .= selectdim(
                A, phi_dim, 1
            )
        end
    end

    return new_grid_vec, new_A
end

function get_interpolator(
        ::Type{<:StructuredGrid},
        A::AbstractArray{T, 4}, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1
    ) where {T}
    As = reinterpret(reshape, SVector{3, T}, A)
    return get_interpolator(StructuredGrid, As, gridr, gridθ, gridϕ, order, bc)
end

function _create_spherical_vector_field_interpolator(interpr, interpθ, interpϕ)
    function get_field(xu)
        r_val, θ_val, ϕ_val = cart2sph(xu)

        Br = interpr(r_val, θ_val, ϕ_val)
        Bθ = interpθ(r_val, θ_val, ϕ_val)
        Bϕ = interpϕ(r_val, θ_val, ϕ_val)

        Bvec = sph_to_cart_vector(Br, Bθ, Bϕ, θ_val, ϕ_val)

        return Bvec
    end
    return get_field
end

function _create_spherical_vector_field_interpolator(itp)
    function get_field(xu)
        r_val, θ_val, ϕ_val = cart2sph(xu)

        B_local = itp(r_val, θ_val, ϕ_val)
        Br, Bθ, Bϕ = B_local

        Bvec = sph_to_cart_vector(Br, Bθ, Bϕ, θ_val, ϕ_val)

        return Bvec
    end
    return get_field
end

function get_interpolator(
        ::Type{<:StructuredGrid},
        A::AbstractArray{T, 3}, gridr, gridθ, gridϕ, order::Int = 1, bc::Int = 1
    ) where {T}
    # Detect if uniform grid (old Spherical) or non-uniform r (old SphericalNonUniformR)
    # We check if gridr is an AbstractRange
    # Base.LogRange is an AbstractRange but not uniform!
    is_uniform_r = gridr isa AbstractRange && !(gridr isa Base.LogRange)

    if is_uniform_r
        itp_unscaled = _get_interp_object(StructuredGrid, A, order, bc)
        itp = scale(itp_unscaled, gridr, gridθ, gridϕ)
    else # Non-uniform R (SphericalNonUniformR behavior)
        bctype = (Flat(), Flat(), Periodic())
        gridϕ, A = _ensure_full_phi(gridϕ, A)
        itp = extrapolate(interpolate!((gridr, gridθ, gridϕ), A, Gridded(Linear())), bctype)
    end

    if T <: SVector
        return _create_spherical_vector_field_interpolator(itp)
    else
        return _create_spherical_scalar_field_interpolator(itp)
    end
end

function _create_spherical_scalar_field_interpolator(interp)
    function get_field(xu)
        r_val, θ_val, ϕ_val = cart2sph(xu[1], xu[2], xu[3])
        return interp(r_val, θ_val, ϕ_val)
    end
    return get_field
end

function _get_interp_object(A, order::Int, bc::Int)
    bspline = _get_bspline(order, bc == 2)

    bctype = if bc == 1
        if eltype(A) <: SVector
            SVector{3, eltype(eltype(A))}(NaN, NaN, NaN)
        else
            eltype(eltype(A))(NaN)
        end
    elseif bc == 2
        Periodic()
    else
        Flat()
    end

    return extrapolate(interpolate(A, bspline), bctype)
end

function _get_interp_object(::Type{<:StructuredGrid}, A, order::Int, bc::Int)
    bspline_r = _get_bspline(order, false)
    bspline_θ = _get_bspline(order, false)
    bspline_ϕ = _get_bspline(order, true)

    itp_type = (bspline_r, bspline_θ, bspline_ϕ)

    bctype = if eltype(A) <: SVector
        SVector{3, eltype(eltype(A))}(NaN, NaN, NaN)
    else
        eltype(A)(NaN)
    end

    return extrapolate(interpolate(A, itp_type), bctype)
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
