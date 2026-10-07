# GPU-compatible Cartesian grid field interpolators.

"""
    GPUGrid3D{T, B, V, A<:AbstractArray{V, 3}} <: AbstractFieldInterpolator

Uniform 3D Cartesian grid interpolator compatible with GPU device execution.
"""
struct GPUGrid3D{T, B, V, A <: AbstractArray{V, 3}} <: AbstractFieldInterpolator
    data::A
    x0::T
    dx::T
    inv_dx::T
    nx::Int32
    xmax::T
    y0::T
    dy::T
    inv_dy::T
    ny::Int32
    ymax::T
    z0::T
    dz::T
    inv_dz::T
    nz::Int32
    zmax::T
    bc::B
end

function GPUGrid3D(
        data::A,
        x0::T, dx::T, inv_dx::T, nx::Integer,
        y0::T, dy::T, inv_dy::T, ny::Integer,
        z0::T, dz::T, inv_dz::T, nz::Integer,
        bc::B
    ) where {T, B, V, A <: AbstractArray{V, 3}}
    xmax = x0 + T(nx - 1) * dx
    ymax = y0 + T(ny - 1) * dy
    zmax = z0 + T(nz - 1) * dz
    return GPUGrid3D(
        data,
        x0, dx, inv_dx, Int32(nx), xmax,
        y0, dy, inv_dy, Int32(ny), ymax,
        z0, dz, inv_dz, Int32(nz), zmax,
        bc
    )
end

Adapt.adapt_structure(to, g::GPUGrid3D) = GPUGrid3D(
    Adapt.adapt(to, g.data),
    g.x0, g.dx, g.inv_dx, g.nx, g.xmax,
    g.y0, g.dy, g.inv_dy, g.ny, g.ymax,
    g.z0, g.dz, g.inv_dz, g.nz, g.zmax,
    Adapt.adapt(to, g.bc),
)

"""
    GPUGrid2D{T, B, V, A<:AbstractArray{V, 2}} <: AbstractFieldInterpolator

Uniform 2D Cartesian grid interpolator compatible with GPU device execution.
"""
struct GPUGrid2D{T, B, V, A <: AbstractArray{V, 2}} <: AbstractFieldInterpolator
    data::A
    x0::T
    dx::T
    inv_dx::T
    nx::Int32
    xmax::T
    y0::T
    dy::T
    inv_dy::T
    ny::Int32
    ymax::T
    bc::B
end

function GPUGrid2D(
        data::A,
        x0::T, dx::T, inv_dx::T, nx::Integer,
        y0::T, dy::T, inv_dy::T, ny::Integer,
        bc::B
    ) where {T, B, V, A <: AbstractArray{V, 2}}
    xmax = x0 + T(nx - 1) * dx
    ymax = y0 + T(ny - 1) * dy
    return GPUGrid2D(
        data,
        x0, dx, inv_dx, Int32(nx), xmax,
        y0, dy, inv_dy, Int32(ny), ymax,
        bc
    )
end

Adapt.adapt_structure(to, g::GPUGrid2D) = GPUGrid2D(
    Adapt.adapt(to, g.data),
    g.x0, g.dx, g.inv_dx, g.nx, g.xmax,
    g.y0, g.dy, g.inv_dy, g.ny, g.ymax,
    Adapt.adapt(to, g.bc),
)

"""
    GPUGrid1D{T, B, V, A<:AbstractArray{V, 1}} <: AbstractFieldInterpolator

Uniform 1D Cartesian grid interpolator compatible with GPU device execution.
"""
struct GPUGrid1D{T, B, V, A <: AbstractArray{V, 1}} <: AbstractFieldInterpolator
    data::A
    x0::T
    dx::T
    inv_dx::T
    nx::Int32
    xmax::T
    dir::Int32
    bc::B
end

function GPUGrid1D(
        data::A,
        x0::T, dx::T, inv_dx::T, nx::Integer,
        dir::Integer,
        bc::B
    ) where {T, B, V, A <: AbstractArray{V, 1}}
    xmax = x0 + T(nx - 1) * dx
    return GPUGrid1D(
        data,
        x0, dx, inv_dx, Int32(nx), xmax,
        Int32(dir),
        bc
    )
end

Adapt.adapt_structure(to, g::GPUGrid1D) = GPUGrid1D(
    Adapt.adapt(to, g.data),
    g.x0, g.dx, g.inv_dx, g.nx, g.xmax,
    g.dir,
    Adapt.adapt(to, g.bc),
)

@inline function _fill_value(bc, ::Type{V}) where {V}
    if hasproperty(bc, :fill_value)
        val = bc.fill_value
        if val isa V
            return val
        elseif val isa Number && V <: SVector
            T = eltype(V)
            return SVector{length(V), T}(ntuple(_ -> T(val), length(V)))
        elseif val isa Number && V <: Number
            return V(val)
        end
    end
    if V <: SVector
        T = eltype(V)
        return SVector{length(V), T}(ntuple(_ -> T(NaN), length(V)))
    elseif V <: Number
        return V(NaN)
    else
        return zero(V)
    end
end

@inline function _fill_value(bc::Tuple, ::Type{V}) where {V}
    return _fill_value(bc[1], V)
end

@inline function _apply_bc_1d(x, x0, xmax, bc)
    isnan(x) && return (x, true)
    if bc isa FillExtrap
        (x < x0 || x > xmax) && return (x, true)
    elseif bc isa ClampExtrap
        x = clamp(x, x0, xmax)
    elseif bc isa WrapExtrap
        span = xmax - x0
        x = mod(x - x0, span) + x0
    end
    return (x, false)
end

@inline function _get_bc_dim(bc::Tuple, dim::Int)
    return bc[dim]
end

@inline function _get_bc_dim(bc, dim::Int)
    return bc
end

@inline _lerp(c0, c1, w) = muladd(w, c1 - c0, c0)

@inline function (g::GPUGrid3D{T, B, V})(x::Real, y::Real, z::Real) where {T, B, V}
    bc_x = _get_bc_dim(g.bc, 1)
    bc_y = _get_bc_dim(g.bc, 2)
    bc_z = _get_bc_dim(g.bc, 3)

    xmax = g.xmax
    ymax = g.ymax
    zmax = g.zmax

    x_adj, out_x = _apply_bc_1d(x, g.x0, xmax, bc_x)
    y_adj, out_y = _apply_bc_1d(y, g.y0, ymax, bc_y)
    z_adj, out_z = _apply_bc_1d(z, g.z0, zmax, bc_z)

    (out_x || out_y || out_z) && return _fill_value(g.bc, V)

    fx = clamp((x_adj - g.x0) * g.inv_dx, zero(T), T(g.nx - Int32(1)))
    fy = clamp((y_adj - g.y0) * g.inv_dy, zero(T), T(g.ny - Int32(1)))
    fz = clamp((z_adj - g.z0) * g.inv_dz, zero(T), T(g.nz - Int32(1)))

    ix0 = min(unsafe_trunc(Int32, fx), g.nx - Int32(2))
    iy0 = min(unsafe_trunc(Int32, fy), g.ny - Int32(2))
    iz0 = min(unsafe_trunc(Int32, fz), g.nz - Int32(2))

    wx = fx - T(ix0)
    wy = fy - T(iy0)
    wz = fz - T(iz0)

    ix = ix0 + Int32(1)
    iy = iy0 + Int32(1)
    iz = iz0 + Int32(1)

    @inbounds begin
        c000 = g.data[ix, iy, iz]
        c100 = g.data[ix + 1, iy, iz]
        c010 = g.data[ix, iy + 1, iz]
        c110 = g.data[ix + 1, iy + 1, iz]
        c001 = g.data[ix, iy, iz + 1]
        c101 = g.data[ix + 1, iy, iz + 1]
        c011 = g.data[ix, iy + 1, iz + 1]
        c111 = g.data[ix + 1, iy + 1, iz + 1]
    end

    c00 = _lerp(c000, c100, wx)
    c10 = _lerp(c010, c110, wx)
    c01 = _lerp(c001, c101, wx)
    c11 = _lerp(c011, c111, wx)

    c0 = _lerp(c00, c10, wy)
    c1 = _lerp(c01, c11, wy)

    return _lerp(c0, c1, wz)
end

@inline (g::GPUGrid3D)(coords::Tuple) = g(coords[1], coords[2], coords[3])
@inline (g::GPUGrid3D)(xu::AbstractVector) = g(xu[1], xu[2], xu[3])
@inline (g::GPUGrid3D)(xu::AbstractVector, t) = g(xu)

@inline function (g::GPUGrid2D{T, B, V})(x::Real, y::Real) where {T, B, V}
    bc_x = _get_bc_dim(g.bc, 1)
    bc_y = _get_bc_dim(g.bc, 2)

    xmax = g.xmax
    ymax = g.ymax

    x_adj, out_x = _apply_bc_1d(x, g.x0, xmax, bc_x)
    y_adj, out_y = _apply_bc_1d(y, g.y0, ymax, bc_y)

    (out_x || out_y) && return _fill_value(g.bc, V)

    fx = clamp((x_adj - g.x0) * g.inv_dx, zero(T), T(g.nx - Int32(1)))
    fy = clamp((y_adj - g.y0) * g.inv_dy, zero(T), T(g.ny - Int32(1)))

    ix0 = min(unsafe_trunc(Int32, fx), g.nx - Int32(2))
    iy0 = min(unsafe_trunc(Int32, fy), g.ny - Int32(2))

    wx = fx - T(ix0)
    wy = fy - T(iy0)

    ix = ix0 + Int32(1)
    iy = iy0 + Int32(1)

    @inbounds begin
        c00 = g.data[ix, iy]
        c10 = g.data[ix + 1, iy]
        c01 = g.data[ix, iy + 1]
        c11 = g.data[ix + 1, iy + 1]
    end

    c0 = _lerp(c00, c10, wx)
    c1 = _lerp(c01, c11, wx)

    return _lerp(c0, c1, wy)
end

@inline (g::GPUGrid2D)(coords::Tuple) = g(coords[1], coords[2])
@inline (g::GPUGrid2D)(xu::AbstractVector) = g(xu[1], xu[2])
@inline (g::GPUGrid2D)(xu::AbstractVector, t) = g(xu)

@inline function (g::GPUGrid1D{T, B, V})(x::Real) where {T, B, V}
    bc_x = _get_bc_dim(g.bc, 1)
    xmax = g.xmax

    x_adj, out_x = _apply_bc_1d(x, g.x0, xmax, bc_x)
    out_x && return _fill_value(g.bc, V)

    fx = clamp((x_adj - g.x0) * g.inv_dx, zero(T), T(g.nx - Int32(1)))
    ix0 = min(unsafe_trunc(Int32, fx), g.nx - Int32(2))
    wx = fx - T(ix0)
    ix = ix0 + Int32(1)

    @inbounds begin
        c0 = g.data[ix]
        c1 = g.data[ix + 1]
    end

    return _lerp(c0, c1, wx)
end

@inline (g::GPUGrid1D)(coords::Tuple) = g(coords[1])
@inline (g::GPUGrid1D)(xu::AbstractVector) = g(xu[g.dir])
@inline (g::GPUGrid1D)(xu::AbstractVector, t) = g(xu)

@inline function _grid_props(g)
    if hasproperty(g, :lo)
        return g.lo, g.h, g.inv_h, g.len
    elseif hasproperty(g, :inner)
        x0 = g.inner[1]
        len = length(g.inner)
        h = g.h isa AbstractArray ? g.h[1] : g.h
        inv_h = g.inv_h isa AbstractArray ? g.inv_h[1] : g.inv_h
        return x0, h, inv_h, len
    elseif hasproperty(g, :x)
        return _grid_props(g.x)
    else
        x0 = first(g)
        len = length(g)
        h = step(g)
        inv_h = inv(h)
        return x0, h, inv_h, len
    end
end

function _to_gpu_grid(itp, backend::Backend; dir = 1)
    if isdefined(itp, :grids) && isdefined(itp, :data)
        N = length(itp.grids)
        if N == 3
            gx, gy, gz = itp.grids
            x0, dx, inv_dx, nx = _grid_props(gx)
            y0, dy, inv_dy, ny = _grid_props(gy)
            z0, dz, inv_dz, nz = _grid_props(gz)
            T = typeof(x0)
            data_gpu = Adapt.adapt(backend, itp.data)
            bc = itp.extraps
            return GPUGrid3D(
                data_gpu,
                T(x0), T(dx), T(inv_dx), Int32(nx),
                T(y0), T(dy), T(inv_dy), Int32(ny),
                T(z0), T(dz), T(inv_dz), Int32(nz),
                bc,
            )
        elseif N == 2
            gx, gy = itp.grids
            x0, dx, inv_dx, nx = _grid_props(gx)
            y0, dy, inv_dy, ny = _grid_props(gy)
            T = typeof(x0)
            data_gpu = Adapt.adapt(backend, itp.data)
            bc = itp.extraps
            return GPUGrid2D(
                data_gpu,
                T(x0), T(dx), T(inv_dx), Int32(nx),
                T(y0), T(dy), T(inv_dy), Int32(ny),
                bc,
            )
        elseif N == 1
            gx = itp.grids[1]
            x0, dx, inv_dx, nx = _grid_props(gx)
            T = typeof(x0)
            data_gpu = Adapt.adapt(backend, itp.data)
            bc = itp.extraps
            return GPUGrid1D(
                data_gpu,
                T(x0), T(dx), T(inv_dx), Int32(nx),
                Int32(dir),
                bc,
            )
        end
    elseif (isdefined(itp, :grid) || isdefined(itp, :x)) &&
            (isdefined(itp, :data) || isdefined(itp, :y))
        g = isdefined(itp, :grid) ? itp.grid : itp.x
        x0, dx, inv_dx, nx = _grid_props(g)
        T = typeof(x0)
        data = isdefined(itp, :data) ? itp.data : itp.y
        data_gpu = Adapt.adapt(backend, data)
        bc = if hasproperty(itp, :extraps)
            itp.extraps
        elseif hasproperty(itp, :extrap)
            itp.extrap
        else
            nothing
        end
        return GPUGrid1D(
            data_gpu,
            T(x0), T(dx), T(inv_dx), Int32(nx),
            Int32(dir),
            bc,
        )
    end
    return Adapt.adapt(backend, itp)
end

"""
    GPUUniformAxis{T, B}

Uniform 1D grid axis for GPU-compatible interpolators.
"""
struct GPUUniformAxis{T, B}
    x0::T
    dx::T
    inv_dx::T
    n::Int32
    xmax::T
    bc::B
end

function GPUUniformAxis(x0::T, dx::T, inv_dx::T, n::Integer, xmax::T, bc::B) where {T, B}
    return GPUUniformAxis{T, B}(x0, dx, inv_dx, Int32(n), xmax, bc)
end

function GPUUniformAxis(x0::T, dx::T, inv_dx::T, n::Integer, bc::B) where {T, B}
    xmax = x0 + T(n - 1) * dx
    return GPUUniformAxis(x0, dx, inv_dx, Int32(n), xmax, bc)
end

Adapt.adapt_structure(to, ax::GPUUniformAxis) = GPUUniformAxis(
    ax.x0, ax.dx, ax.inv_dx, ax.n, ax.xmax, Adapt.adapt(to, ax.bc),
)

"""
    GPUNonUniformAxis{T, V<:AbstractVector{T}, B}

Non-uniform 1D grid axis for GPU-compatible interpolators.
"""
struct GPUNonUniformAxis{T, V <: AbstractVector{T}, B}
    coords::V
    x0::T
    xmax::T
    n::Int32
    bc::B
end

function GPUNonUniformAxis(
        coords::V, x0::T, xmax::T, n::Integer, bc::B
    ) where {T, V <: AbstractVector{T}, B}
    return GPUNonUniformAxis{T, V, B}(coords, x0, xmax, Int32(n), bc)
end

function GPUNonUniformAxis(coords::V, bc::B) where {T, V <: AbstractVector{T}, B}
    x0 = coords[1]
    xmax = coords[end]
    n = Int32(length(coords))
    return GPUNonUniformAxis{T, V, B}(coords, x0, xmax, n, bc)
end

Adapt.adapt_structure(to, ax::GPUNonUniformAxis) = GPUNonUniformAxis(
    Adapt.adapt(to, ax.coords), ax.x0, ax.xmax, ax.n, Adapt.adapt(to, ax.bc),
)

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
    v_adj, out = _apply_bc_1d(v, ax.x0, ax.xmax, ax.bc)
    out && return (Int32(1), zero(T), true)
    fv = clamp((v_adj - ax.x0) * ax.inv_dx, zero(T), T(ax.n - Int32(1)))
    i0 = min(unsafe_trunc(Int32, fv), ax.n - Int32(2))
    w = fv - T(i0)
    return (i0 + Int32(1), w, false)
end

@inline function _axis_eval(ax::GPUNonUniformAxis{T}, val::Real) where {T}
    v = T(val)
    v_adj, out = _apply_bc_1d(v, ax.x0, ax.xmax, ax.bc)
    out && return (Int32(1), zero(T), true)
    i, w = _search_bracket(ax.coords, v_adj, ax.n)
    return (i, w, false)
end

"""
    GPUSphericalGrid{T, Ar, Aθ, Aϕ, V, A, B} <: AbstractFieldInterpolator

Spherical grid interpolator compatible with GPU device execution.
Supports both uniform and non-uniform grids in r, θ, and ϕ.
"""
struct GPUSphericalGrid{
    T, Ar, Aθ, Aϕ, V, A <: AbstractArray{V, 3}, B
} <: AbstractFieldInterpolator
    data::A
    axis_r::Ar
    axis_θ::Aθ
    axis_ϕ::Aϕ
    bc::B
end

function GPUSphericalGrid(
        data::A,
        axis_r::Ar,
        axis_θ::Aθ,
        axis_ϕ::Aϕ,
        bc::B
    ) where {Ar, Aθ, Aϕ, V, A <: AbstractArray{V, 3}, B}
    T = typeof(axis_r.x0)
    return GPUSphericalGrid{T, Ar, Aθ, Aϕ, V, A, B}(data, axis_r, axis_θ, axis_ϕ, bc)
end

Adapt.adapt_structure(to, g::GPUSphericalGrid) = GPUSphericalGrid(
    Adapt.adapt(to, g.data),
    Adapt.adapt(to, g.axis_r),
    Adapt.adapt(to, g.axis_θ),
    Adapt.adapt(to, g.axis_ϕ),
    Adapt.adapt(to, g.bc),
)

@inline _transform_spherical_result(res::SVector{3}, θ, ϕ) =
    sph2cartvec(res[1], res[2], res[3], θ, ϕ)
@inline _transform_spherical_result(res::Number, θ, ϕ) = res
@inline _transform_spherical_result(res, θ, ϕ) =
    length(res) > 1 ? sph2cartvec(res[1], res[2], res[3], θ, ϕ) : res

@inline function (g::GPUSphericalGrid{T, Ar, Aθ, Aϕ, V})(
        x::Real, y::Real, z::Real
    ) where {T, Ar, Aθ, Aϕ, V}
    rθϕ = cart2sph(x, y, z)
    r = rθϕ[1]
    θ = rθϕ[2]
    ϕ = rθϕ[3]

    ir, wr, out_r = _axis_eval(g.axis_r, r)
    iθ, wθ, out_θ = _axis_eval(g.axis_θ, θ)
    iϕ, wϕ, out_ϕ = _axis_eval(g.axis_ϕ, ϕ)

    (out_r || out_θ || out_ϕ) && return _fill_value(g.bc, V)

    @inbounds begin
        c000 = g.data[ir, iθ, iϕ]
        c100 = g.data[ir + 1, iθ, iϕ]
        c010 = g.data[ir, iθ + 1, iϕ]
        c110 = g.data[ir + 1, iθ + 1, iϕ]
        c001 = g.data[ir, iθ, iϕ + 1]
        c101 = g.data[ir + 1, iθ, iϕ + 1]
        c011 = g.data[ir, iθ + 1, iϕ + 1]
        c111 = g.data[ir + 1, iθ + 1, iϕ + 1]
    end

    c00 = _lerp(c000, c100, wr)
    c10 = _lerp(c010, c110, wr)
    c01 = _lerp(c001, c101, wr)
    c11 = _lerp(c011, c111, wr)

    c0 = _lerp(c00, c10, wθ)
    c1 = _lerp(c01, c11, wθ)

    res = _lerp(c0, c1, wϕ)
    return _transform_spherical_result(res, θ, ϕ)
end

@inline (g::GPUSphericalGrid)(coords::Tuple) = g(coords[1], coords[2], coords[3])
@inline (g::GPUSphericalGrid)(coords::Tuple, t) = g(coords[1], coords[2], coords[3])
@inline (g::GPUSphericalGrid)(xu::AbstractVector) = g(xu[1], xu[2], xu[3])
@inline (g::GPUSphericalGrid)(xu::AbstractVector, t) = g(xu[1], xu[2], xu[3])
@inline (g::GPUSphericalGrid)(x::Real, y::Real, z::Real, t) = g(x, y, z)

function _to_gpu_axis(g, backend::Backend, bc)
    if hasproperty(g, :lo) && g.h isa Number
        T = typeof(g.lo)
        x0 = T(g.lo)
        dx = T(g.h)
        inv_dx = T(g.inv_h)
        n = Int32(g.len)
        xmax = x0 + T(n - 1) * dx
        return GPUUniformAxis(x0, dx, inv_dx, n, xmax, bc)
    elseif hasproperty(g, :inner)
        coords = g.inner
        T = eltype(coords)
        coords_gpu = Adapt.adapt(backend, coords)
        n = Int32(length(coords))
        return GPUNonUniformAxis(coords_gpu, T(coords[1]), T(coords[end]), n, bc)
    elseif g isa AbstractRange && !(g isa Base.LogRange)
        T = eltype(g)
        x0 = T(first(g))
        dx = T(step(g))
        inv_dx = inv(dx)
        n = Int32(length(g))
        xmax = x0 + T(n - 1) * dx
        return GPUUniformAxis(x0, dx, inv_dx, n, xmax, bc)
    elseif g isa AbstractVector
        T = eltype(g)
        coords_gpu = Adapt.adapt(backend, collect(g))
        n = Int32(length(g))
        return GPUNonUniformAxis(coords_gpu, T(first(g)), T(last(g)), n, bc)
    else
        throw(ArgumentError("Unsupported grid axis type: $(typeof(g))"))
    end
end

function _to_gpu_spherical_grid(itp, backend::Backend)
    if isdefined(itp, :grids) && isdefined(itp, :data)
        gr, gθ, gϕ = itp.grids
        bc = itp.extraps
        bcr = _get_bc_dim(bc, 1)
        bcθ = _get_bc_dim(bc, 2)
        bcϕ = _get_bc_dim(bc, 3)

        axis_r = _to_gpu_axis(gr, backend, bcr)
        axis_θ = _to_gpu_axis(gθ, backend, bcθ)
        axis_ϕ = _to_gpu_axis(gϕ, backend, bcϕ)

        data_gpu = Adapt.adapt(backend, itp.data)
        return GPUSphericalGrid(data_gpu, axis_r, axis_θ, axis_ϕ, bc)
    else
        throw(ArgumentError("Expected interpolant with grids and data fields."))
    end
end

