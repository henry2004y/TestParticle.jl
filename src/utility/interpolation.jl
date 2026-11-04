# Field interpolations.

"""
    AbstractFieldInterpolator

Abstract type for all field interpolators.
"""
abstract type AbstractFieldInterpolator <: Function end

# The device conversion needs the grids, boundary conditions and order an interpolant was
# built from. Mutable so that interpolators hold it by reference: solvers copy the field
# parameters by value into several structs, and an inline spec makes every copy larger.
mutable struct InterpSpec{G <: Tuple, B <: Tuple}
    const grids::G
    const bcs::B
    const order::Int
end

Adapt.adapt_structure(to, s::InterpSpec) =
    InterpSpec(Adapt.adapt(to, s.grids), s.bcs, s.order)

# Vector grids are copied so that a caller mutating theirs cannot desynchronize the spec
# from the interpolant.
_own_grid(g::AbstractRange) = g
_own_grid(g::AbstractVector) = collect(g)

"""
    FieldInterpolator{T, S}

A callable struct that wraps a 3D interpolation object.
"""
struct FieldInterpolator{T, S <: InterpSpec} <: AbstractFieldInterpolator
    itp::T
    spec::S
end

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

Adapt.adapt_structure(to, fi::FieldInterpolator) =
    FieldInterpolator(Adapt.adapt(to, fi.itp), Adapt.adapt(to, fi.spec))

"""
    FieldInterpolator2D{T, S}

A callable struct that wraps a 2D interpolation object.
"""
struct FieldInterpolator2D{T, S <: InterpSpec} <: AbstractFieldInterpolator
    itp::T
    spec::S
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

Adapt.adapt_structure(to, fi::FieldInterpolator2D) =
    FieldInterpolator2D(Adapt.adapt(to, fi.itp), Adapt.adapt(to, fi.spec))

"""
    FieldInterpolator1D{T, S}

A callable struct that wraps a 1D interpolation object.
"""
struct FieldInterpolator1D{T, S <: InterpSpec} <: AbstractFieldInterpolator
    itp::T
    spec::S
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

Adapt.adapt_structure(to, fi::FieldInterpolator1D) =
    FieldInterpolator1D(Adapt.adapt(to, fi.itp), Adapt.adapt(to, fi.spec), fi.dir)

function jacobian(f, x)
    return ForwardDiff.jacobian(f, x)
end

"""
    SphericalFieldInterpolator{T, S}

A callable struct for spherical grid interpolation (scalar or combined vector).
"""
struct SphericalFieldInterpolator{T, S <: InterpSpec} <: AbstractFieldInterpolator
    itp::T
    spec::S
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

Adapt.adapt_structure(to, fi::SphericalFieldInterpolator) =
    SphericalFieldInterpolator(Adapt.adapt(to, fi.itp), Adapt.adapt(to, fi.spec))

function _cardinal_interp(grids::Tuple, A, coeffs, bc, extrap, store)
    itp = cardinal_interp(grids, A; coeffs, bc, extrap, store)
    itp.data isa Array && return itp

    @warn "FastInterpolations cannot evaluate ND cardinal interpolation with " *
        "aliased non-Array data; falling back to owned storage." maxlog = 1
    return cardinal_interp(
        grids, A; coeffs, bc, extrap, store = StorePolicy()
    )
end

@inline function _cardinal_interp(grid, A, coeffs, bc, extrap, store)
    return cardinal_interp(grid, A; coeffs, bc, extrap, store)
end

function _fastinterp(
        grids, A, order, bc_or_extrap, coeffs = OnTheFly();
        store::StorePolicy = StorePolicy()
    )
    axes = grids isa Tuple ? grids : (grids,)
    for g in axes
        length(g) >= 2 && all(isfinite, g) && all(d -> d > zero(d), diff(g)) ||
            throw(ArgumentError("Grid nodes must be finite and strictly increasing."))
    end
    policies = bc_or_extrap isa Tuple ? bc_or_extrap : map(Returns(bc_or_extrap), axes)
    length(policies) == length(axes) ||
        throw(ArgumentError("Provide one boundary policy per grid axis."))
    policies = map(policies) do policy
        if policy isa FillExtrap && policy.fill_value isa Number
            T = eltype(A)
            value = T <: SVector ? T(ntuple(_ -> policy.fill_value, length(T))) :
                T(policy.fill_value)
            return FillExtrap(value)
        end
        return policy
    end
    bcs = map(p -> p isa AbstractExtrap ? NoBC() : p, policies)
    extrap = map(p -> p isa AbstractExtrap ? p : NoExtrap(), policies)
    bc_val = grids isa Tuple ? bcs : only(bcs)
    extrap_val = grids isa Tuple ? extrap : only(extrap)

    itp = if order == 1
        linear_interp(grids, A; bc = bc_val, extrap = extrap_val, store)
    elseif order == 3
        _cardinal_interp(grids, A, coeffs, bc_val, extrap_val, store)
    elseif order == 0
        constant_interp(grids, A; bc = bc_val, extrap = extrap_val, store)
    else
        throw(ArgumentError("Interpolation order $order is not supported. Supported orders are 0, 1, and 3."))
    end
    axes = grids isa Tuple ? grids : (grids,)
    spec = InterpSpec(map(_own_grid, axes), policies, order)
    return itp, spec
end

# Not public API: FastInterpolations has no accessor for the stored node values.
_itp_data(itp) = hasproperty(itp, :data) ? itp.data : itp.y

_device_bc(bc::Union{FillExtrap, ClampExtrap, WrapExtrap}) = bc
# An inclusive period repeats the first node at the end, i.e. wraps over the grid span.
_device_bc(::PeriodicBC{:inclusive}) = WrapExtrap()
_device_bc(bc) = throw(
    ArgumentError(
        "Boundary condition $bc has no device implementation. Use FillExtrap, " *
            "ClampExtrap, WrapExtrap or an inclusive PeriodicBC."
    )
)

function _device_bcs(spec::InterpSpec)
    spec.order == 1 || throw(
        ArgumentError(
            "Device grids interpolate linearly, so an interpolator of order " *
                "$(spec.order) cannot be moved to a device. Build it with order = 1."
        )
    )
    return map(_device_bc, spec.bcs)
end

@inline function _apply_extrap(x, x0::T, xmax::T, ::FillExtrap) where {T}
    # Rounding in a coordinate transform, e.g. `cart2sph`, may place a point
    # marginally outside the grid it was built from, so widen by a few ulps.
    tol = 4 * eps(T) * max(abs(x0), abs(xmax))
    return (x, x < x0 - tol || x > xmax + tol)
end

@inline _apply_extrap(x, x0, xmax, ::ClampExtrap) = (clamp(x, x0, xmax), false)

@inline function _apply_extrap(x, x0, xmax, ::WrapExtrap)
    return (mod(x - x0, xmax - x0) + x0, false)
end

_fill_value(::Tuple{}, ::Type{V}) where {V} = _nan_fill(V)
_fill_value(bcs::Tuple, ::Type{V}) where {V} = first(bcs) isa FillExtrap ?
    _fill_value(first(bcs), V) : _fill_value(Base.tail(bcs), V)

function _fill_value(bc::FillExtrap, ::Type{V}) where {V}
    val = bc.fill_value
    val isa V && return val
    val isa Number || return _nan_fill(V)
    return V <: SVector ? V(ntuple(_ -> val, length(V))) : V(val)
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

"""
    build_interpolator(A, grids::Tuple; bc, order=1, gridtype=CartesianGrid, kwargs...)
    build_interpolator([gridtype,] A, grids..., order=1, bc=FillExtrap(NaN); kwargs...)

Build a callable field for [`prepare`](@ref). Vector data has shape
`(3, length.(grids)...)` or contains `SVector{3}` values. Orders are 0, 1, and 3.
The tuple API requires `bc`: one policy for all axes or a tuple of per-axis policies.

`PeriodicBC(endpoint=:inclusive)` repeats the endpoint; `:exclusive` omits it
(uniform-grid period `NΔx`). `WrapExtrap()` wraps over the grid span without imposing
spline periodicity. `FillExtrap(NaN)` fills outside; `ClampExtrap()` extends constantly.
Field lookup neither wraps particle state nor terminates trajectories.

`gridtype=StructuredGrid` uses spherical `(r, θ, ϕ)` grids and `(Br, Bθ, Bϕ)` components;
queries and outputs are Cartesian. Supply fill/clamp for radius and colatitude and
`PeriodicBC` with period `2π` for longitude. Its positional default is
`(FillExtrap(NaN), FillExtrap(NaN), PeriodicBC(check = false))`.

`coeffs` selects cubic coefficient storage. `store=StorePolicy(copy=false)` may alias
inputs, which must remain unchanged. Device conversion requires order 1 and fill,
clamp, wrap, or inclusive periodic policies; exclusive periodic fields are CPU-only.
"""
function build_interpolator end

function build_interpolator(
        A::AbstractArray, grids::Tuple; bc, order::Int = 1,
        gridtype = CartesianGrid, kwargs...
    )
    1 <= length(grids) <= 3 || throw(ArgumentError("Provide 1, 2, or 3 grid axes."))
    return build_interpolator(gridtype, A, grids..., order, bc; kwargs...)
end

@inline build_interpolator(A::AbstractArray, grid1, args...; kwargs...) =
    build_interpolator(CartesianGrid, A, grid1, args...; kwargs...)

function build_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc = FillExtrap(NaN);
        coeffs = OnTheFly(), store::StorePolicy = StorePolicy()
    ) where {T}
    @assert size(A, 1) == 3 "Incompatible 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(
        CartesianGrid, As, gridx, gridy, gridz, order, bc; coeffs, store
    )
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc = FillExtrap(NaN);
        coeffs = OnTheFly(), store::StorePolicy = StorePolicy()
    ) where {T}
    _check_interpolation_consistency(A, (gridx, gridy, gridz), order)
    itp, spec = _fastinterp((gridx, gridy, gridz), A, order, bc, coeffs; store)
    return FieldInterpolator(itp, spec)
end

function build_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 4},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc = FillExtrap(NaN);
        store::StorePolicy = StorePolicy()
    ) where {T}
    @assert size(A, 1) == 3 "Incompatible 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(
        RectilinearGrid, As, gridx, gridy, gridz, order, bc; store
    )
end

function build_interpolator(
        ::Type{<:RectilinearGrid}, A::AbstractArray{T, 3},
        gridx::AbstractVector, gridy::AbstractVector, gridz::AbstractVector,
        order::Int = 1, bc = FillExtrap(NaN);
        store::StorePolicy = StorePolicy()
    ) where {T}
    if order != 1
        throw(
            ArgumentError(
                "RectilinearGrid (CartesianNonUniform) only supports " *
                    "order=1 (Linear) interpolation."
            )
        )
    end

    itp, spec = _fastinterp((gridx, gridy, gridz), A, order, bc; store)
    return FieldInterpolator(itp, spec)
end

function build_interpolator(
        ::Type{<:StructuredGrid}, A::AbstractArray{T, 4},
        gridr, gridθ, gridϕ, order::Int = 1,
        bc = (FillExtrap(NaN), FillExtrap(NaN), PeriodicBC(check = false));
        coeffs = OnTheFly(), store::StorePolicy = StorePolicy()
    ) where {T}
    @assert size(A, 1) == 3 "Incompatible 3D force field and grid!"
    As = reinterpret(reshape, SVector{3, T}, A)
    return build_interpolator(
        StructuredGrid, As, gridr, gridθ, gridϕ, order, bc; coeffs, store
    )
end

function build_interpolator(
        ::Type{<:StructuredGrid}, A::AbstractArray{T, 3},
        gridr, gridθ, gridϕ, order::Int = 1,
        bc = (FillExtrap(NaN), FillExtrap(NaN), PeriodicBC(check = false));
        coeffs = OnTheFly(), store::StorePolicy = StorePolicy()
    ) where {T}
    r_min = minimum(gridr)
    θ_min, θ_max = extrema(gridθ)
    ϕ_min, ϕ_max = extrema(gridϕ)

    @assert r_min >= 0 "r must be non-negative."
    @assert θ_min >= 0 && θ_max <= Float64(π) + 4eps(float(θ_max)) "θ must be within [0, π]."
    @assert ϕ_min >= 0 && ϕ_max <= 2π + 4eps(float(ϕ_max)) "ϕ must be within [0, 2π]."

    bc isa Tuple && length(bc) == 3 || throw(
        ArgumentError(
            "Spherical fields require (radial, colatitude, longitude) boundary policies."
        )
    )
    all(p -> p isa Union{FillExtrap, ClampExtrap}, bc[1:2]) ||
        throw(ArgumentError("Radius and colatitude support fill or clamp policies."))
    longitude = bc[3]
    longitude isa PeriodicBC || throw(
        ArgumentError(
            "Longitude requires PeriodicBC with an explicit endpoint convention."
        )
    )
    span = last(gridϕ) - first(gridϕ)
    period = if longitude isa PeriodicBC{:inclusive}
        span
    elseif isnothing(longitude.period)
        steps = diff(gridϕ)
        all(d -> isapprox(d, first(steps)), steps) || throw(
            ArgumentError(
                "Nonuniform exclusive longitude requires an explicit period."
            )
        )
        span + first(steps)
    else
        longitude.period
    end
    isapprox(period, oftype(float(period), 2π)) ||
        throw(ArgumentError("Longitude period must be 2π."))
    _check_interpolation_consistency(A, (gridr, gridθ, gridϕ), order)
    itp, spec = _fastinterp((gridr, gridθ, gridϕ), A, order, bc, coeffs; store)
    return SphericalFieldInterpolator(itp, spec)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A,
        gridx::AbstractVector, gridy::AbstractVector,
        order::Int = 1, bc = FillExtrap(NaN);
        coeffs = OnTheFly(), store::StorePolicy = StorePolicy()
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 2 "Incompatible 2D force field and grid! " *
            "Expected 2D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 3 "Incompatible 2D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    _check_interpolation_consistency(As, (gridx, gridy), order)
    itp, spec = _fastinterp((gridx, gridy), As, order, bc, coeffs; store)
    return FieldInterpolator2D(itp, spec)
end

function build_interpolator(
        ::Type{<:CartesianGrid}, A, gridx::AbstractVector,
        order::Int = 1, bc = FillExtrap(NaN);
        dir = 1, coeffs = OnTheFly(), store::StorePolicy = StorePolicy()
    )
    if eltype(A) <: SVector
        @assert ndims(A) == 1 "Incompatible 1D force field and grid! " *
            "Expected 1D array of SVectors."
        As = A
    else
        @assert size(A, 1) == 3 && ndims(A) == 2 "Incompatible 1D force field and grid!"
        As = reinterpret(reshape, SVector{3, eltype(A)}, A)
    end

    dir in 1:3 || throw(ArgumentError("dir must be 1, 2, or 3."))
    _check_interpolation_consistency(As, (gridx,), order)
    itp, spec = _fastinterp(gridx, As, order, bc, coeffs; store)
    return FieldInterpolator1D(itp, spec, dir)
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
