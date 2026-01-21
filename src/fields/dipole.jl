# Magnetic dipole field.

@kwdef struct DipoleField{M} <: AbstractField{false}
    BMoment::M = BMoment_Earth
end

(M::DipoleField)(xu) = dipole(xu[1:3], M.BMoment)

"""
    getB_dipole(xu, BMoment = BMoment_Earth)

Return dipole magnetic field in [T].
"""
@inline function getB_dipole(xu; BMoment = BMoment_Earth)
    return dipole(xu[1:3], BMoment)
end

"""
Calculates the magnetic field from a dipole with magnetic moment `M` at `r`.
"""
@inline function dipole(rIn::AbstractVector, M::SVector{3})
    x, y, z = rIn
    r = sqrt(x^2 + y^2 + z^2)
    Coef = μ₀ / (4 * π * r^5)

    m11 = 3 * x^2 - r^2
    m12 = 3 * x * y
    m13 = 3 * x * z
    m21 = 3 * y * x
    m22 = 3 * y^2 - r^2
    m23 = 3 * y * z
    m31 = 3 * z * x
    m32 = 3 * z * y
    m33 = 3 * z^2 - r^2

    mat = SMatrix{3, 3}(m11, m21, m31, m12, m22, m32, m13, m23, m33)
    return mat * M * Coef
end

"""
    dipole_fieldline(ϕ, L=2.5, nP=100)

Creates `nP` points on one field line of the magnetic field from a dipole. In a centered
dipole magnetic field model, the path along a given L shell can be described as r = L*cos²λ,
where r is the radial distance (in planetary radii) to a point on the line,
λ is its co-latitude, and L is the L-shell of interest.
"""
function dipole_fieldline(ϕ, L = 2.5, nP::Int = 100)
    xyz = [sph2cart(L * sin(θ)^2, θ, ϕ) for θ in range(0, stop = π, length = nP)]
    x = getindex.(xyz, 1)
    y = getindex.(xyz, 2)
    z = getindex.(xyz, 3)
    return (x, y, z)
end
