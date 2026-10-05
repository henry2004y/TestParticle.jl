# 3D Morton Z-order curve particle ordering for GPU cache optimization.
# Reference: Morton (1966), "A computer Oriented Geodetic Data Base; and a New
# Technique in File Sequencing".

@inline function _part1by2(n::UInt32)
    n &= 0x000003ff
    n = (n | (n << 16)) & 0x030000ff
    n = (n | (n << 8)) & 0x0300f00f
    n = (n | (n << 4)) & 0x030c30c3
    n = (n | (n << 2)) & 0x09249249
    return n
end

"""
    morton3D(ix::Integer, iy::Integer, iz::Integer) -> UInt32

Compute the 3D Morton code (Z-order) for 10-bit integer coordinates `(ix, iy, iz)`.
Input coordinates are clamped to `[0, 1023]`.
"""
@inline function morton3D(ix::Integer, iy::Integer, iz::Integer)
    return (_part1by2(UInt32(clamp(ix, 0, 1023))) << 2) |
        (_part1by2(UInt32(clamp(iy, 0, 1023))) << 1) |
        _part1by2(UInt32(clamp(iz, 0, 1023)))
end

"""
    morton_sort_particles(xv_init::AbstractMatrix{T}) -> Vector{Int}

Compute a permutation vector that sorts particles along the 3D Morton curve
based on their initial positions `(x, y, z) = xv_init[:, 1:3]`.
"""
function morton_sort_particles(xv_init::AbstractMatrix{T}) where {T}
    N = size(xv_init, 1)
    xmin, xmax = extrema(@view xv_init[:, 1])
    ymin, ymax = extrema(@view xv_init[:, 2])
    zmin, zmax = extrema(@view xv_init[:, 3])

    sx = xmax > xmin ? T(1023) / (xmax - xmin) : zero(T)
    sy = ymax > ymin ? T(1023) / (ymax - ymin) : zero(T)
    sz = zmax > zmin ? T(1023) / (zmax - zmin) : zero(T)

    codes = Vector{UInt32}(undef, N)
    @inbounds for i in 1:N
        x_val, y_val, z_val = xv_init[i, 1], xv_init[i, 2], xv_init[i, 3]
        ix = isnan(x_val) ? 0 : clamp(floor(Int, (x_val - xmin) * sx), 0, 1023)
        iy = isnan(y_val) ? 0 : clamp(floor(Int, (y_val - ymin) * sy), 0, 1023)
        iz = isnan(z_val) ? 0 : clamp(floor(Int, (z_val - zmin) * sz), 0, 1023)
        codes[i] = morton3D(ix, iy, iz)
    end
    return sortperm(codes)
end
