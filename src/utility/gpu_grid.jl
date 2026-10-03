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
    y0::T
    dy::T
    inv_dy::T
    ny::Int32
    z0::T
    dz::T
    inv_dz::T
    nz::Int32
    bc::B
end

Adapt.adapt_structure(to, g::GPUGrid3D) = GPUGrid3D(
    Adapt.adapt(to, g.data),
    g.x0, g.dx, g.inv_dx, g.nx,
    g.y0, g.dy, g.inv_dy, g.ny,
    g.z0, g.dz, g.inv_dz, g.nz,
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
    y0::T
    dy::T
    inv_dy::T
    ny::Int32
    bc::B
end

Adapt.adapt_structure(to, g::GPUGrid2D) = GPUGrid2D(
    Adapt.adapt(to, g.data),
    g.x0, g.dx, g.inv_dx, g.nx,
    g.y0, g.dy, g.inv_dy, g.ny,
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
    dir::Int32
    bc::B
end

Adapt.adapt_structure(to, g::GPUGrid1D) = GPUGrid1D(
    Adapt.adapt(to, g.data),
    g.x0, g.dx, g.inv_dx, g.nx,
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

@inline function (g::GPUGrid3D{T, B, V})(x::Real, y::Real, z::Real) where {T, B, V}
    bc_x = _get_bc_dim(g.bc, 1)
    bc_y = _get_bc_dim(g.bc, 2)
    bc_z = _get_bc_dim(g.bc, 3)

    xmax = g.x0 + T(g.nx - Int32(1)) * g.dx
    ymax = g.y0 + T(g.ny - Int32(1)) * g.dy
    zmax = g.z0 + T(g.nz - Int32(1)) * g.dz

    x_adj, out_x = _apply_bc_1d(x, g.x0, xmax, bc_x)
    y_adj, out_y = _apply_bc_1d(y, g.y0, ymax, bc_y)
    z_adj, out_z = _apply_bc_1d(z, g.z0, zmax, bc_z)

    (out_x || out_y || out_z) && return _fill_value(g.bc, V)

    fx = clamp((x_adj - g.x0) * g.inv_dx, zero(T), T(g.nx - Int32(1)))
    fy = clamp((y_adj - g.y0) * g.inv_dy, zero(T), T(g.ny - Int32(1)))
    fz = clamp((z_adj - g.z0) * g.inv_dz, zero(T), T(g.nz - Int32(1)))

    ix = min(floor(Int32, fx) + Int32(1), g.nx - Int32(1))
    iy = min(floor(Int32, fy) + Int32(1), g.ny - Int32(1))
    iz = min(floor(Int32, fz) + Int32(1), g.nz - Int32(1))

    wx = fx - T(ix - Int32(1))
    wy = fy - T(iy - Int32(1))
    wz = fz - T(iz - Int32(1))

    c000 = g.data[ix, iy, iz]
    c100 = g.data[ix + 1, iy, iz]
    c010 = g.data[ix, iy + 1, iz]
    c110 = g.data[ix + 1, iy + 1, iz]
    c001 = g.data[ix, iy, iz + 1]
    c101 = g.data[ix + 1, iy, iz + 1]
    c011 = g.data[ix, iy + 1, iz + 1]
    c111 = g.data[ix + 1, iy + 1, iz + 1]

    c00 = c000 * (one(T) - wx) + c100 * wx
    c10 = c010 * (one(T) - wx) + c110 * wx
    c01 = c001 * (one(T) - wx) + c101 * wx
    c11 = c011 * (one(T) - wx) + c111 * wx

    c0 = c00 * (one(T) - wy) + c10 * wy
    c1 = c01 * (one(T) - wy) + c11 * wy

    return c0 * (one(T) - wz) + c1 * wz
end

@inline (g::GPUGrid3D)(coords::Tuple) = g(coords[1], coords[2], coords[3])
@inline (g::GPUGrid3D)(xu::AbstractVector) = g(xu[1], xu[2], xu[3])
@inline (g::GPUGrid3D)(xu::AbstractVector, t) = g(xu)

@inline function (g::GPUGrid2D{T, B, V})(x::Real, y::Real) where {T, B, V}
    bc_x = _get_bc_dim(g.bc, 1)
    bc_y = _get_bc_dim(g.bc, 2)

    xmax = g.x0 + T(g.nx - Int32(1)) * g.dx
    ymax = g.y0 + T(g.ny - Int32(1)) * g.dy

    x_adj, out_x = _apply_bc_1d(x, g.x0, xmax, bc_x)
    y_adj, out_y = _apply_bc_1d(y, g.y0, ymax, bc_y)

    (out_x || out_y) && return _fill_value(g.bc, V)

    fx = clamp((x_adj - g.x0) * g.inv_dx, zero(T), T(g.nx - Int32(1)))
    fy = clamp((y_adj - g.y0) * g.inv_dy, zero(T), T(g.ny - Int32(1)))

    ix = min(floor(Int32, fx) + Int32(1), g.nx - Int32(1))
    iy = min(floor(Int32, fy) + Int32(1), g.ny - Int32(1))

    wx = fx - T(ix - Int32(1))
    wy = fy - T(iy - Int32(1))

    c00 = g.data[ix, iy]
    c10 = g.data[ix + 1, iy]
    c01 = g.data[ix, iy + 1]
    c11 = g.data[ix + 1, iy + 1]

    c0 = c00 * (one(T) - wx) + c10 * wx
    c1 = c01 * (one(T) - wx) + c11 * wx

    return c0 * (one(T) - wy) + c1 * wy
end

@inline (g::GPUGrid2D)(coords::Tuple) = g(coords[1], coords[2])
@inline (g::GPUGrid2D)(xu::AbstractVector) = g(xu[1], xu[2])
@inline (g::GPUGrid2D)(xu::AbstractVector, t) = g(xu)

@inline function (g::GPUGrid1D{T, B, V})(x::Real) where {T, B, V}
    bc_x = _get_bc_dim(g.bc, 1)
    xmax = g.x0 + T(g.nx - Int32(1)) * g.dx

    x_adj, out_x = _apply_bc_1d(x, g.x0, xmax, bc_x)
    out_x && return _fill_value(g.bc, V)

    fx = clamp((x_adj - g.x0) * g.inv_dx, zero(T), T(g.nx - Int32(1)))
    ix = min(floor(Int32, fx) + Int32(1), g.nx - Int32(1))
    wx = fx - T(ix - Int32(1))

    c0 = g.data[ix]
    c1 = g.data[ix + 1]

    return c0 * (one(T) - wx) + c1 * wx
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
    elseif (isdefined(itp, :grid) || isdefined(itp, :x)) && isdefined(itp, :data)
        g = isdefined(itp, :grid) ? itp.grid : itp.x
        x0, dx, inv_dx, nx = _grid_props(g)
        T = typeof(x0)
        data_gpu = Adapt.adapt(backend, itp.data)
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
