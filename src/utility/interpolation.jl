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
`FillExtrap` values must agree across axes.

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

# Examples
```jldoctest
julia> grid = (-1.0:1.0:2.0, -1.0:1.0:2.0, -1.0:1.0:2.0);

julia> B_arr = zeros(3, 4, 4, 4);

julia> B_arr[3, :, :, :] .= 1.0e-8;

julia> B = build_interpolator(B_arr, grid; bc = TestParticle.FillExtrap(NaN));

julia> p = prepare(; B, species = Proton);

julia> p.B(SA[0.5, 0.5, 0.5], 0.0) ≈ SA[0.0, 0.0, 1.0e-8]
true
```
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

# Wrapped in a mutable struct with const fields so atomic slot swaps exchange a single
# pointer even when `F` is an isbits callable or tuple.
mutable struct _TimeSlab{T, F}
    const idx::Int
    const idx2::Int
    const t1::T
    const t2::T
    const inv_dt::T
    const f1::F
    const f2::F
end

mutable struct _AtomicSlab{T, F}
    @atomic slab::Union{Nothing, _TimeSlab{T, F}}
end

mutable struct _LRUNode{F}
    const idx::Int
    const field::F
    prev::Union{Nothing, _LRUNode{F}}
    next::Union{Nothing, _LRUNode{F}}
end

mutable struct _LRUCache{F}
    const capacity::Int
    const map::Dict{Int, _LRUNode{F}}
    head::Union{Nothing, _LRUNode{F}}
    tail::Union{Nothing, _LRUNode{F}}
end

function _detach!(cache::_LRUCache{F}, node::_LRUNode{F}) where {F}
    p = node.prev
    n = node.next
    if p === nothing
        cache.head = n
    else
        p.next = n
    end
    if n === nothing
        cache.tail = p
    else
        n.prev = p
    end
    node.prev = nothing
    node.next = nothing
    return node
end

function _push_tail!(cache::_LRUCache{F}, node::_LRUNode{F}) where {F}
    t = cache.tail
    node.prev = t
    node.next = nothing
    if t === nothing
        cache.head = node
    else
        t.next = node
    end
    cache.tail = node
    return node
end

function _lru_get!(cache::_LRUCache{F}, idx::Int) where {F}
    node = get(cache.map, idx, nothing)
    node === nothing && return nothing
    if node !== cache.tail
        _detach!(cache, node)
        _push_tail!(cache, node)
    end
    return node
end

function _invalidate_slots!(
        slots::Vector{_AtomicSlab{T, F}}, evicted_idx::Int
    ) where {T, F}
    @inbounds for s in slots
        slab = @atomic :acquire s.slab
        if slab !== nothing && (slab.idx == evicted_idx || slab.idx2 == evicted_idx)
            @atomic :release s.slab = nothing
        end
    end
    return nothing
end

function _lru_put!(
        cache::_LRUCache{F}, idx::Int, field, slots::Vector{_AtomicSlab{T, F}}
    ) where {T, F}
    existing = get(cache.map, idx, nothing)
    if existing !== nothing
        if existing !== cache.tail
            _detach!(cache, existing)
            _push_tail!(cache, existing)
        end
        return existing
    end
    if length(cache.map) >= cache.capacity
        victim = cache.head
        if victim !== nothing
            _detach!(cache, victim)
            delete!(cache.map, victim.idx)
            _invalidate_slots!(slots, victim.idx)
        end
    end
    node = _LRUNode{F}(idx, field, nothing, nothing)
    cache.map[idx] = node
    return _push_tail!(cache, node)
end

_own_time_grid(times::AbstractRange) = float(times)
_own_time_grid(times::AbstractVector) = collect(float(eltype(times)), times)

function _validate_times(t_grid::AbstractVector{T}) where {T <: AbstractFloat}
    n = length(t_grid)
    n >= 2 || throw(ArgumentError("Time grid must contain at least 2 points."))
    @inbounds begin
        all(isfinite, t_grid) ||
            throw(ArgumentError("Time grid points must be finite."))
        for i in 1:(n - 1)
            t_grid[i + 1] > t_grid[i] ||
                throw(ArgumentError("Time grid points must be strictly increasing."))
        end
    end
    return t_grid
end

function _resolve_time_end(t_grid::AbstractVector{T}, bc) where {T}
    t_first = first(t_grid)
    t_last = last(t_grid)
    if bc isa Union{ClampExtrap, FillExtrap, WrapExtrap}
        return t_last
    elseif bc isa PeriodicBC{:inclusive}
        isnothing(bc.period) || throw(
            ArgumentError("Inclusive PeriodicBC does not accept an explicit period.")
        )
        return t_last
    elseif bc isa PeriodicBC{:exclusive}
        span = t_last - t_first
        period = if isnothing(bc.period)
            if t_grid isa AbstractRange
                span + T(step(t_grid))
            else
                dt0 = @inbounds t_grid[2] - t_grid[1]
                @inbounds for i in 2:(length(t_grid) - 1)
                    isapprox(t_grid[i + 1] - t_grid[i], dt0) || throw(
                        ArgumentError(
                            "Nonuniform exclusive PeriodicBC requires an explicit period."
                        )
                    )
                end
                span + dt0
            end
        else
            bc.period isa Real && isfinite(bc.period) ||
                throw(ArgumentError("PeriodicBC period must be finite."))
            p = T(bc.period)
            p > span || throw(
                ArgumentError("Exclusive PeriodicBC period must exceed time span.")
            )
            p
        end
        return t_first + period
    else
        throw(
            ArgumentError(
                "Unsupported time boundary policy $bc. Use ClampExtrap, " *
                    "FillExtrap, WrapExtrap, or PeriodicBC."
            )
        )
    end
end

"""
    LazyTimeInterpolator(times, loader; cache_size = ..., bc = ClampExtrap())

Callable interpolator for time-dependent fields with lazy loading and linear interpolation
between consecutive snapshots on `times`.

# Arguments
- `times::AbstractVector`: Strictly increasing finite real time nodes (`length >= 2`).
- `loader`: Callable `i -> field` loading snapshot `i` (or a `Tuple`/`NamedTuple` of fields
  for shared multi-field caching via unpacking, e.g. `E, B = LazyTimeInterpolator(...)`).
- `cache_size::Int`: Maximum LRU snapshot capacity (`>= 2`, default
  `max(2, 2 * Threads.nthreads())`).
- `bc`: Time boundary policy (`ClampExtrap()`, `FillExtrap(val)`, `WrapExtrap()`, or
  `PeriodicBC()`).

# Examples
```jldoctest
julia> times = 0.0:1.0:2.0;

julia> loader(i) = x -> SA[x[1] * times[i], x[2], x[3]];

julia> itp = TestParticle.LazyTimeInterpolator(times, loader);

julia> itp(SA[2.0, 0.0, 0.0], 0.5) ≈ SA[1.0, 0.0, 0.0]
true
```
"""
struct LazyTimeInterpolator{T, F, L, V <: AbstractVector{T}, B} <: Function
    times::V
    t_first::T
    t_last::T
    t_end::T
    loader::L
    bc::B
    slots::Vector{_AtomicSlab{T, F}}
    slot_mask::Int
    cache::_LRUCache{F}
    lock::ReentrantLock
end

function LazyTimeInterpolator(
        times::AbstractVector, loader;
        cache_size::Int = max(2, 2 * Threads.nthreads()),
        bc = ClampExtrap()
    )
    eltype(times) <: Real || throw(ArgumentError("Time grid points must be real."))
    cache_size >= 2 || throw(ArgumentError("cache_size must be at least 2."))

    t_grid = _validate_times(_own_time_grid(times))
    t_first = first(t_grid)
    t_last = last(t_grid)
    t_end = _resolve_time_end(t_grid, bc)

    f1 = loader(1)
    if f1 isa Union{Tuple, NamedTuple}
        isempty(f1) &&
            throw(ArgumentError("Multi-field loader cannot return empty fields."))
    end

    return _LazyTimeInterpolator(t_grid, t_first, t_last, t_end, loader, bc, cache_size, f1)
end

function _LazyTimeInterpolator(
        t_grid::AbstractVector{T}, t_first::T, t_last::T, t_end::T,
        loader, bc, cache_size::Int, f1::F
    ) where {T, F}
    nslots = nextpow(2, max(16, 2 * Threads.maxthreadid()))
    slots = [_AtomicSlab{T, F}(nothing) for _ in 1:nslots]
    cache = _LRUCache{F}(cache_size, Dict{Int, _LRUNode{F}}(), nothing, nothing)
    _lru_put!(cache, 1, f1, slots)
    return LazyTimeInterpolator{T, F, typeof(loader), typeof(t_grid), typeof(bc)}(
        t_grid, t_first, t_last, t_end, loader, bc,
        slots, nslots - 1, cache, ReentrantLock()
    )
end

struct _WholeField end
struct _ComponentSel{K} end
struct _OpValue end
struct _OpJacobian end
struct _OpDerivT end

# Component view returned when indexing or unpacking a multi-field LazyTimeInterpolator so
# E and B can be passed separately to prepare(; E, B) while sharing one snapshot cache.
struct _LazyFieldComponent{K, I <: LazyTimeInterpolator} <: Function
    parent::I
end

function Base.length(::LazyTimeInterpolator{T, F}) where {T, F <: Union{Tuple, NamedTuple}}
    return fieldcount(F)
end

@inline function Base.getindex(
        itp::LazyTimeInterpolator{T, F}, k::Int
    ) where {T, F <: Union{Tuple, NamedTuple}}
    1 <= k <= fieldcount(F) || throw(BoundsError(itp, k))
    return _LazyFieldComponent{k, typeof(itp)}(itp)
end

@inline function Base.iterate(
        itp::LazyTimeInterpolator{T, F}, state::Int = 1
    ) where {T, F <: Union{Tuple, NamedTuple}}
    state > fieldcount(F) && return nothing
    return (itp[state], state + 1)
end

@inline function Base.indexed_iterate(
        itp::LazyTimeInterpolator{T, F}, i::Int, state::Int = 1
    ) where {T, F <: Union{Tuple, NamedTuple}}
    return (itp[i], i + 1)
end

Base.keys(::LazyTimeInterpolator{T, F}) where {T, F <: NamedTuple} = fieldnames(F)

@inline function Base.getindex(
        itp::LazyTimeInterpolator{T, F}, sym::Symbol
    ) where {T, F <: NamedTuple}
    k = Base.fieldindex(F, sym, false)
    k == 0 && throw(KeyError(sym))
    return _LazyFieldComponent{k, typeof(itp)}(itp)
end

@inline function Base.getproperty(
        itp::LazyTimeInterpolator{T, F}, sym::Symbol
    ) where {T, F <: NamedTuple}
    hasfield(typeof(itp), sym) && return getfield(itp, sym)
    return itp[sym]
end

function Base.propertynames(itp::LazyTimeInterpolator{T, F}) where {T, F <: NamedTuple}
    return (fieldnames(typeof(itp))..., fieldnames(F)...)
end

@inline _select(::_WholeField, f) = f
@inline _select(::_ComponentSel{K}, f) where {K} = getfield(f, K)

@inline _call_field(f, x) = f(x)
@inline _call_field(fs::Tuple, x) = map(f -> f(x), fs)
@inline _call_field(fs::NamedTuple{K}, x) where {K} =
    NamedTuple{K}(_call_field(Tuple(fs), x))

@inline _call_jacobian(f, x) = jacobian(f, x)
@inline _call_jacobian(fs::Tuple, x) = map(f -> jacobian(f, x), fs)
@inline _call_jacobian(fs::NamedTuple{K}, x) where {K} =
    NamedTuple{K}(_call_jacobian(Tuple(fs), x))

@inline _lincomb(y1, y2, w) = (one(w) - w) * y1 + w * y2
@inline _lincomb(y1::Tuple, y2::Tuple, w) = map((a, b) -> _lincomb(a, b, w), y1, y2)
@inline _lincomb(y1::NamedTuple{K}, y2::NamedTuple{K}, w) where {K} =
    NamedTuple{K}(_lincomb(Tuple(y1), Tuple(y2), w))

@inline _scale_diff(y2, y1, inv_dt) = (y2 - y1) * inv_dt
@inline _scale_diff(y2::Tuple, y1::Tuple, inv_dt) =
    map((b, a) -> _scale_diff(b, a, inv_dt), y2, y1)
@inline _scale_diff(y2::NamedTuple{K}, y1::NamedTuple{K}, inv_dt) where {K} =
    NamedTuple{K}(_scale_diff(Tuple(y2), Tuple(y1), inv_dt))

@inline _zero_like(v) = zero(v)
@inline _zero_like(v::Tuple) = map(_zero_like, v)
@inline _zero_like(v::NamedTuple{K}) where {K} = NamedTuple{K}(_zero_like(Tuple(v)))

@inline _has_nan(v::Number) = isnan(v)
@inline _has_nan(v) = any(_has_nan, v)

@inline function _fill_value_for(v::V, fill_val) where {V <: StaticArray}
    if fill_val isa V
        return fill_val
    elseif fill_val isa AbstractVector && length(fill_val) == length(V)
        return V(fill_val)
    elseif fill_val isa Number
        E = eltype(V)
        return V(ntuple(_ -> E(fill_val), Val(length(V))))
    else
        E = eltype(V)
        return V(ntuple(_ -> E(NaN), Val(length(V))))
    end
end

@inline function _fill_value_for(v::V, fill_val) where {V <: Number}
    return fill_val isa Number ? V(fill_val) : V(NaN)
end

@inline function _fill_value_for(v::AbstractArray, fill_val)
    E = eltype(v)
    fill_el = fill_val isa Number ? E(fill_val) : E(NaN)
    out = similar(v)
    fill!(out, fill_el)
    return out
end

@inline function _fill_value_for(v::Tuple, fill_val)
    if fill_val isa Tuple && length(fill_val) == length(v)
        return map(_fill_value_for, v, fill_val)
    else
        return map(x -> _fill_value_for(x, fill_val), v)
    end
end

@inline function _fill_value_for(v::NamedTuple{K}, fill_val) where {K}
    return NamedTuple{K}(_fill_value_for(Tuple(v), fill_val))
end

@inline function _eval_slab(::_OpValue, sel, slab::_TimeSlab, x, t)
    w = (t - slab.t1) * slab.inv_dt
    f1 = _select(sel, slab.f1)
    f2 = _select(sel, slab.f2)
    iszero(w) && return _call_field(f1, x)
    isone(w) && return _call_field(f2, x)
    return _lincomb(_call_field(f1, x), _call_field(f2, x), w)
end

@inline function _eval_slab(::_OpJacobian, sel, slab::_TimeSlab, x, t)
    w = (t - slab.t1) * slab.inv_dt
    f1 = _select(sel, slab.f1)
    f2 = _select(sel, slab.f2)
    iszero(w) && return _call_jacobian(f1, x)
    isone(w) && return _call_jacobian(f2, x)
    return _lincomb(_call_jacobian(f1, x), _call_jacobian(f2, x), w)
end

@inline function _eval_slab(::_OpDerivT, sel, slab::_TimeSlab, x, t)
    f1 = _select(sel, slab.f1)
    f2 = _select(sel, slab.f2)
    return _scale_diff(_call_field(f2, x), _call_field(f1, x), slab.inv_dt)
end

@inline _eval_single(::_OpValue, sel, field, x) = _call_field(_select(sel, field), x)
@inline _eval_single(::_OpJacobian, sel, field, x) = _call_jacobian(_select(sel, field), x)
@inline _eval_single(::_OpDerivT, sel, field, x) =
    _zero_like(_call_field(_select(sel, field), x))

@inline function _eval_fill_oob(::_OpValue, sel, field, x, fill_val)
    return _fill_value_for(_call_field(_select(sel, field), x), fill_val)
end

@inline function _eval_fill_oob(::_OpJacobian, sel, field, x, fill_val)
    sample = _call_jacobian(_select(sel, field), x)
    return _has_nan(fill_val) ? _fill_value_for(sample, NaN) : _zero_like(sample)
end

@inline function _eval_fill_oob(::_OpDerivT, sel, field, x, fill_val)
    sample = _call_field(_select(sel, field), x)
    return _has_nan(fill_val) ? _fill_value_for(sample, NaN) : _zero_like(sample)
end

@inline function _current_slot(itp::LazyTimeInterpolator)
    tid = Threads.threadid()
    return @inbounds itp.slots[((tid - 1) & itp.slot_mask) + 1]
end

@inline function _in_slab(
        slab::_TimeSlab{T}, tv::Real, t_end::T, ::Union{ClampExtrap, FillExtrap}
    ) where {T}
    return slab.t1 <= tv && (tv < slab.t2 || (tv == t_end && slab.t2 == t_end))
end

@inline function _in_slab(slab::_TimeSlab, tv::Real, t_end, bc)
    return slab.t1 <= tv && tv < slab.t2
end

@inline function _find_slab_index(times::AbstractRange, tv::Real)
    n = length(times)
    u = (tv - first(times)) / step(times)
    i = if u <= 0
        1
    elseif u >= n - 1
        n - 1
    else
        clamp(unsafe_trunc(Int, floor(u)) + 1, 1, n - 1)
    end
    @inbounds if tv < times[i] && i > 1
        i -= 1
    elseif tv >= times[i + 1] && i < n - 1
        i += 1
    end
    return i
end

@inline function _find_slab_index(times::AbstractVector, tv::Real)
    return clamp(searchsortedlast(times, tv), 1, length(times) - 1)
end

function _fetch_slab!(
        itp::LazyTimeInterpolator{T, F}, slot::_AtomicSlab{T, F},
        idx::Int, idx2::Int, t1::T, t2::T
    ) where {T, F}
    return lock(itp.lock) do
        cur = @atomic :acquire slot.slab
        if cur !== nothing && cur.idx == idx && cur.idx2 == idx2 &&
                cur.t1 == t1 && cur.t2 == t2
            return cur
        end
        # Promote cached endpoints before loading missing ones so cache_size = 2 never
        # evicts a sibling endpoint that is needed for the same slab.
        n1 = _lru_get!(itp.cache, idx)
        n2 = idx2 == idx ? n1 : _lru_get!(itp.cache, idx2)
        if n1 === nothing
            f1 = itp.loader(idx)
            n1 = _lru_put!(itp.cache, idx, f1, itp.slots)
        end
        if n2 === nothing
            f2 = itp.loader(idx2)
            n2 = _lru_put!(itp.cache, idx2, f2, itp.slots)
        end
        slab = _TimeSlab{T, F}(idx, idx2, t1, t2, inv(t2 - t1), n1.field, n2.field)
        @atomic :release slot.slab = slab
        return slab
    end
end

@inline function _get_slab!(
        itp::LazyTimeInterpolator, slot::_AtomicSlab, idx::Int, idx2::Int
    )
    t1 = @inbounds itp.times[idx]
    t2 = @inbounds itp.times[idx2]
    return _fetch_slab!(itp, slot, idx, idx2, t1, t2)
end

@inline function _get_single_field!(
        itp::LazyTimeInterpolator, slab::Union{Nothing, _TimeSlab}, idx::Int
    )
    if slab !== nothing
        slab.idx == idx && return slab.f1
        slab.idx2 == idx && return slab.f2
    end
    return _get_single_field_locked!(itp, idx)
end

function _get_single_field_locked!(itp::LazyTimeInterpolator, idx::Int)
    return lock(itp.lock) do
        node = _lru_get!(itp.cache, idx)
        if node === nothing
            field = itp.loader(idx)
            node = _lru_put!(itp.cache, idx, field, itp.slots)
        end
        return node.field
    end
end

function _eval_lazy_slow(
        op, sel, itp::LazyTimeInterpolator{T, F, L, V, <:ClampExtrap},
        slot::_AtomicSlab, cached_slab::Union{Nothing, _TimeSlab}, x, t, tv::Real
    ) where {T, F, L, V}
    if tv < itp.t_first
        return _eval_single(op, sel, _get_single_field!(itp, cached_slab, 1), x)
    elseif tv > itp.t_last
        n = length(itp.times)
        return _eval_single(op, sel, _get_single_field!(itp, cached_slab, n), x)
    end
    idx = tv == itp.t_last ? length(itp.times) - 1 : _find_slab_index(itp.times, tv)
    slab = _get_slab!(itp, slot, idx, idx + 1)
    return _eval_slab(op, sel, slab, x, t)
end

function _eval_lazy_slow(
        op, sel, itp::LazyTimeInterpolator{T, F, L, V, <:FillExtrap},
        slot::_AtomicSlab, cached_slab::Union{Nothing, _TimeSlab}, x, t, tv::Real
    ) where {T, F, L, V}
    tol = 4 * eps(T) * max(abs(itp.t_first), abs(itp.t_last))
    if tv < itp.t_first - tol || tv > itp.t_last + tol
        field1 = _get_single_field!(itp, cached_slab, 1)
        return _eval_fill_oob(op, sel, field1, x, itp.bc.fill_value)
    end
    t_q = if tv < itp.t_first
        zero(t) + itp.t_first
    elseif tv > itp.t_last
        zero(t) + itp.t_last
    else
        t
    end
    tv_q = ForwardDiff.value(t_q)
    idx = tv_q == itp.t_last ? length(itp.times) - 1 : _find_slab_index(itp.times, tv_q)
    slab = _get_slab!(itp, slot, idx, idx + 1)
    return _eval_slab(op, sel, slab, x, t_q)
end

function _eval_lazy_slow(
        op, sel,
        itp::LazyTimeInterpolator{
            T, F, L, V, <:Union{WrapExtrap, PeriodicBC{:inclusive}, PeriodicBC{:exclusive}},
        },
        slot::_AtomicSlab, cached_slab::Union{Nothing, _TimeSlab}, x, t, tv::Real
    ) where {T, F, L, V}
    period = itp.t_end - itp.t_first
    if itp.t_first <= tv < itp.t_end
        t_q = t
        tv_q = tv
    else
        t_q = mod(t - itp.t_first, period) + itp.t_first
        tv_q = ForwardDiff.value(t_q)
        if tv_q >= itp.t_end
            t_q -= period
            tv_q = ForwardDiff.value(t_q)
        end
        if tv_q < itp.t_first
            t_q = zero(t) + itp.t_first
            tv_q = itp.t_first
        end
    end
    if cached_slab !== nothing && cached_slab.t1 <= tv_q < cached_slab.t2
        return _eval_slab(op, sel, cached_slab, x, t_q)
    end
    slab = if itp.bc isa PeriodicBC{:exclusive} && tv_q >= itp.t_last
        _fetch_slab!(itp, slot, length(itp.times), 1, itp.t_last, itp.t_end)
    else
        idx = _find_slab_index(itp.times, tv_q)
        _get_slab!(itp, slot, idx, idx + 1)
    end
    return _eval_slab(op, sel, slab, x, t_q)
end

@inline function _eval_lazy(op, sel, itp::LazyTimeInterpolator, x, t)
    tv = ForwardDiff.value(t)
    slot = _current_slot(itp)
    slab = @atomic :acquire slot.slab
    if isnan(tv)
        return _eval_fill_oob(op, sel, _get_single_field!(itp, slab, 1), x, NaN)
    end
    if slab !== nothing && _in_slab(slab, tv, itp.t_end, itp.bc)
        return _eval_slab(op, sel, slab, x, t)
    end
    return _eval_lazy_slow(op, sel, itp, slot, slab, x, t, tv)
end

@inline (itp::LazyTimeInterpolator)(x, t) = _eval_lazy(_OpValue(), _WholeField(), itp, x, t)
@inline function (c::_LazyFieldComponent{K})(x, t) where {K}
    return _eval_lazy(_OpValue(), _ComponentSel{K}(), c.parent, x, t)
end

@inline function jacobian(itp::LazyTimeInterpolator, x, t)
    return _eval_lazy(_OpJacobian(), _WholeField(), itp, x, t)
end
@inline function jacobian(c::_LazyFieldComponent{K}, x, t) where {K}
    return _eval_lazy(_OpJacobian(), _ComponentSel{K}(), c.parent, x, t)
end

@inline function derivative_t(itp::LazyTimeInterpolator, x, t)
    return _eval_lazy(_OpDerivT(), _WholeField(), itp, x, t)
end
@inline function derivative_t(c::_LazyFieldComponent{K}, x, t) where {K}
    return _eval_lazy(_OpDerivT(), _ComponentSel{K}(), c.parent, x, t)
end
