# Collection of utility functions and commonly used constants.

include("constants.jl")
include("zero.jl")

"""
    sph2cart(r, θ, ϕ)

Convert from spherical to Cartesian coordinates vector.
"""
function sph2cart(r, θ, ϕ)
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    return SVector{3}(r * sinθ * cosϕ, r * sinθ * sinϕ, r * cosθ)
end

@inline @inbounds sph2cart(x) = sph2cart(x[1], x[2], x[3])

"""
    cart2sph(x, y, z)

Convert from Cartesian to spherical coordinates vector.
"""
function cart2sph(x, y, z)
    r = hypot(x, y, z)
    if r == 0
        return SVector{3, eltype(r)}(0, 0, 0)
    end
    θ = acos(z / r)
    ϕ = atan(y, x) |> mod2pi
    return SVector{3}(r, θ, ϕ)
end

@inline @inbounds cart2sph(x) = cart2sph(x[1], x[2], x[3])

"""
    sph_to_cart_vector(vr, vθ, vϕ, θ, ϕ)

Convert a vector from spherical to Cartesian coordinates. θ and ϕ are defined in radians.
"""
function sph_to_cart_vector(vr, vθ, vϕ, θ, ϕ)
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    vx = sinθ * cosϕ * vr + cosθ * cosϕ * vθ - sinϕ * vϕ
    vy = sinθ * sinϕ * vr + cosθ * sinϕ * vθ + cosϕ * vϕ
    vz = cosθ * vr - sinθ * vθ
    return SVector{3}(vx, vy, vz)
end

"""
Return uniform range from 2D/3D CartesianGrid.
"""
function makegrid(grid::CartesianGrid)
    gridmin = coords(minimum(grid))
    gridmax = coords(maximum(grid))
    Δx = spacing(grid)
    dim = paramdim(grid)

    gridx = range(gridmin.x.val, gridmax.x.val, step = Δx[1].val)
    gridy = range(gridmin.y.val, gridmax.y.val, step = Δx[2].val)
    if dim == 3
        gridz = range(gridmin.z.val, gridmax.z.val, step = Δx[3].val)
        return gridx, gridy, gridz
    elseif dim == 2
        return gridx, gridy
    elseif dim == 1
        return (gridx,)
    end
end

"""
Return ranges from 2D/3D Meshes.jl RectilinearGrid.
"""
function makegrid(grid::RectilinearGrid)
    if paramdim(grid) == 3
        return grid.xyz[1], grid.xyz[2], grid.xyz[3]
    elseif paramdim(grid) == 2
        return grid.xyz[1], grid.xyz[2]
    end
end

"""
Return cell center coordinates from 2D/3D CartesianGrid.
"""
function get_cell_centers(grid::CartesianGrid)
    grid_coords = makegrid(grid)
    return map(
        r -> range(first(r) + step(r) / 2, step = step(r), length = length(r) - 1),
        grid_coords
    )
end

"""
Return cell center coordinates from 2D/3D RectilinearGrid.
"""
function get_cell_centers(grid::RectilinearGrid)
    grid_coords = makegrid(grid)
    if grid_coords === nothing
        return nothing
    end
    return map(c -> (c[1:(end - 1)] .+ c[2:end]) ./ 2, grid_coords)
end

"""
    get_perp_vector(b::AbstractVector)

Obtain two unit vectors `e1` and `e2` such that `(e1, e2, b)` form a right-handed orthonormal
system.
"""
function get_perp_vector(b)
    T = eltype(b)
    b̂ = normalize(b)
    v = abs(b̂[3]) < 0.9 ? SA[zero(T), zero(T), one(T)] : SA[zero(T), one(T), zero(T)]
    e1 = normalize(v × b̂)
    e2 = b̂ × e1
    return e1, e2
end

"""
    get_mean_magnitude(B)

Calculate the Root Mean Square (RMS) magnitude of a vector field `B`.
It is assumed that the first dimension of `B` represents the vector components.
This function is compatible with any spatial dimension.
"""
@inline get_mean_magnitude(B) = norm(B) / sqrt(length(B) / size(B, 1))

"""
    get_gyrofrequency(B=5e-9; q=qᵢ, m=mᵢ)

Return the gyrofrequency [rad/s].

# Arguments

  - `B`: Magnetic field magnitude [T]. Default is 5 nT.
  - `q`: Charge [C]. Default is proton charge.
  - `m`: Mass [kg]. Default is proton mass.
"""
get_gyrofrequency(B = 5.0e-9; q = qᵢ, m = mᵢ) = ω = q * B / m

@inline _get_exb_drift(E, B, Bmag) = (E × B) / Bmag^2

@inline _get_v_perp(v, b̂) = v - (v ⋅ b̂) * b̂

function _get_v_gamma(u, func_name)
    is_relativistic = func_name in (
        :trace_relativistic!,
        :trace_relativistic,
        :trace_relativistic_normalized!,
        :trace_relativistic_normalized,
    )

    if is_relativistic
        is_normalized = func_name in (
            :trace_relativistic_normalized!,
            :trace_relativistic_normalized,
        )
        v_real = is_normalized ? get_relativistic_v(u; c = 1) : get_relativistic_v(u)
        v_mag = norm(v_real)
        γ = is_normalized ? 1 / sqrt(1 - v_mag^2) : 1 / sqrt(1 - (v_mag / c)^2)
    else
        v_real = u
        γ = 1.0
    end
    return v_real, γ
end

"""
    get_gyroradius(V, B; q=qᵢ, m=mᵢ)

Return the gyroradius [m].

# Arguments

  - `V`: Velocity magnitude [m/s] (usually perpendicular to the magnetic field).
  - `B`: Magnetic field magnitude [T].
  - `q`: Charge [C]. Default is proton charge.
  - `m`: Mass [kg]. Default is proton mass.
"""
function get_gyroradius(V::Number, B::Number; q = qᵢ, m = mᵢ)
    ω = get_gyrofrequency(B; q, m)
    return r = V / ω
end

function get_gyroradius(v, B, E = SA[0.0, 0.0, 0.0]; q = qᵢ, m = mᵢ)
    Bmag = norm(B)
    if Bmag == 0
        return Inf
    end

    v_E = _get_exb_drift(E, B, Bmag) # ExB drift
    v_gyro = v - v_E

    # Calculate perpendicular component
    b̂ = B / Bmag
    v_gyro_perp = _get_v_perp(v_gyro, b̂)
    V_gyro_perp = norm(v_gyro_perp)

    q2m = q / m
    return V_gyro_perp / abs(q2m * Bmag)
end

"""
    get_gyroradius(sol::AbstractODESolution, t)

Return the gyroradius [m] from the solution `sol` at time `t`.
The interpolated magnetic field function is obtained from `sol.prob.p`.
"""
function get_gyroradius(sol::AbstractODESolution, t)
    # Interpolate state at time t
    xu = sol(t)
    x = get_x(xu)
    u = get_v(xu)

    # Extract parameters
    p = sol.prob.p
    q2m = p[1]
    Efunc = p[3]
    Bfunc = p[4]

    # Calculate fields at position x
    E = Efunc(x, t)
    B = Bfunc(x, t)
    Bmag = norm(B)
    if Bmag == 0
        return Inf
    end

    # Calculate ExB drift
    v_E = _get_exb_drift(E, B, Bmag)

    func_name = Symbol(sol.prob.f.f)
    v_real, γ = _get_v_gamma(u, func_name)

    v_gyro = v_real - v_E # Subtract drift velocity

    # Calculate perpendicular component of v_gyro
    b̂ = B / Bmag
    v_gyro_perp = _get_v_perp(v_gyro, b̂)
    V_gyro_perp = norm(v_gyro_perp)

    # Calculate gyroradius
    # r = γ * m * v_perp / |q * B| = γ * m * v_perp / |q2m * m * B| = γ * v_perp / |q2m * B|
    return γ * V_gyro_perp / abs(q2m * Bmag)
end

"""
    get_gyroperiod(B=5e-9; q=qᵢ, m=mᵢ)

Return the gyroperiod [s].

# Arguments

  - `B`: Magnetic field magnitude [T]. Default is 5 nT.
  - `q`: Charge [C]. Default is proton charge.
  - `m`: Mass [kg]. Default is proton mass.
"""
function get_gyroperiod(B = 5.0e-9; q = qᵢ, m = mᵢ)
    ω = get_gyrofrequency(B; q, m)
    return 2π / ω
end

"""
Return velocity from relativistic γv in `sol`.
"""
function get_velocity(sol)
    v = Array{eltype(sol.u[1]), 2}(undef, 3, length(sol))
    inv_c2 = 1 / c^2
    @inbounds for is in axes(v, 2)
        γv = @view sol[4:6, is]
        γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
        γ = √(1 + γ²v² * inv_c2)
        inv_γ = inv(γ)
        for i in axes(v, 1)
            v[i, is] = γv[i] * inv_γ
        end
    end

    return v
end

"""
Return the energy [eV] from relativistic `sol`.
"""
function get_energy(sol::AbstractODESolution; m = mᵢ, q = qᵢ)
    e = Vector{eltype(sol.u[1])}(undef, length(sol))
    inv_c2 = 1 / c^2
    mc2_q = m * c^2 / abs(q)
    @inbounds for i in eachindex(e)
        γv = @view sol[4:6, i]
        γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
        γ = √(1 + γ²v² * inv_c2)
        e[i] = (γ - 1) * mc2_q
    end

    return e
end

"""
Calculate the energy [eV] of a relativistic particle from γv in [m/s].
"""
function get_energy(γv; m = mᵢ, q = qᵢ)
    γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
    γ = √(1 + γ²v² / c^2)

    return (γ - 1) * m * c^2 / abs(q)
end

"""
Return velocity magnitude from energy in [eV].
"""
energy2velocity(Ek; m = mᵢ, q = qᵢ) = c * sqrt(1 - 1 / (1 + Ek * abs(q) / (m * c^2))^2)

"""
    sample_unit_sphere()

Sample a unit vector on a sphere uniformly.
"""
function sample_unit_sphere()
    ϕ = 2π * rand()
    cosθ = 2 * rand() - 1
    sinθ = sqrt(1 - cosθ^2)
    x = sinθ * cos(ϕ)
    y = sinθ * sin(ϕ)
    z = cosθ

    return SVector{3}(x, y, z)
end

"""
    generate_sphere(nθ=64, nϕ=64, radius=1.0)

Generate a sphere with `radius` and matching grid points `nθ` and `nϕ`.
Returns a tuple of `(x, y, z)` matrices.
"""
function generate_sphere(nθ = 64, nϕ = 64, radius = 1.0)
    ϕ = LinRange(0, 2π, nϕ)
    θ = LinRange(0, π, nθ)
    x = Matrix{Float64}(undef, nϕ, nθ)
    y = Matrix{Float64}(undef, nϕ, nθ)
    z = Matrix{Float64}(undef, nϕ, nθ)
    @inbounds for j in axes(x, 2)
        sinθ, cosθ = sincos(θ[j])
        for i in axes(x, 1)
            sinϕ, cosϕ = sincos(ϕ[i])
            x[i, j] = radius * cosϕ * sinθ
            y[i, j] = radius * sinϕ * sinθ
            z[i, j] = radius * cosθ
        end
    end
    return x, y, z
end

"""
    sample_maxwellian(Tn, m; offset=0.0, u0=SA[0.0, 0.0, 0.0])
    sample_maxwellian(Tn, species::String; offset=0.0, u0=SA[0.0, 0.0, 0.0])

Sample a velocity [m/s] from a Maxwellian distribution with temperature `Tn` [K], mass `m`
[kg], energy offset `offset` [J], and bulk velocity `u0` [m/s].
This function requires `VelocityDistributionFunctions.jl` to be loaded.
"""
function sample_maxwellian end

@inline function _get_curvature(B, Bmag, b̂, JB)
    # ∇|B| = (J_B' * b̂)
    ∇B = JB' * b̂
    # Curvature vector κ = (b̂ ⋅ ∇) b̂
    return (JB * b̂ + b̂ * (-∇B ⋅ b̂)) / Bmag, ∇B
end

@inline function _get_B_jacobian(x, t, Bfunc)
    result = DiffResults.JacobianResult(x)
    result = ForwardDiff.jacobian!(result, r -> Bfunc(r, t), x)
    JB = DiffResults.jacobian(result)
    B = DiffResults.value(result)
    return B, JB
end

"""
    get_magnetic_properties(x, t, Bfunc)

Calculate magnetic field properties at position `x` and time `t`.
Returns tuple `(B, ∇B, κ, b̂, Bmag)`:
- `B`: Magnetic field vector
- `∇B`: Gradient of magnetic field magnitude
- `κ`: Curvature vector
- `b̂`: Unit magnetic field vector
- `Bmag`: Magnitude of B
"""
@inline function get_magnetic_properties(x, t, Bfunc)
    B, JB = _get_B_jacobian(x, t, Bfunc)

    Bmag = norm(B)
    if Bmag == 0
        vzero = zero(x)
        return vzero, vzero, vzero, vzero, Bmag
    end
    b̂ = B / Bmag

    # Share curvature calculation logic
    κ, ∇B = _get_curvature(B, Bmag, b̂, JB)

    return B, ∇B, κ, b̂, Bmag
end

"""
    get_curvature(x, t, Bfunc)

Calculate the curvature vector `κ` of the magnetic field at position `x` and time `t`.
"""
@inline function get_curvature(x, t, Bfunc)
    B, JB = _get_B_jacobian(x, t, Bfunc)

    # Curvature vector κ = (b̂ ⋅ ∇) b̂
    κ = (JB * b̂ + b̂ * (-∇B ⋅ b̂)) / Bmag

    return B, ∇B, κ, b̂, Bmag
end

"""
    get_curvature_radius(x, t, Bfunc)

Calculate the radius of curvature of the magnetic field at position `x` and time `t`.
Returns `Inf` if the field is zero or the field lines are straight.
"""
@inline function get_curvature_radius(x, t, Bfunc)
    _, _, κ, _, _ = get_magnetic_properties(x, t, Bfunc)

    k_mag = norm(κ)
    return iszero(k_mag) ? Inf : inv(k_mag)
end

"""
    get_adiabaticity(r, Bfunc, q, m, μ, t=0.0)
    get_adiabaticity(r, Bfunc, μ, t=0.0; species = Proton)

Calculate the adiabaticity parameter `ϵ = ρ / Rc` at position `r` and time `t`.
`ρ` is the gyroradius and `Rc` is the radius of curvature of the magnetic field.
"""
@inline function get_adiabaticity(r, Bfunc, q, m, μ, t = 0.0)
    Bmag = Bfunc(r, t) |> norm
    iszero(Bmag) && return Inf

    ρ = sqrt(2 * μ * m / Bmag) / abs(q) # Gyroradius

    _, _, κ, _, _ = get_magnetic_properties(r, t, Bfunc)

    k_mag = norm(κ)
    invRc = iszero(k_mag) ? zero(k_mag) : k_mag

    return ρ * invRc
end

function get_adiabaticity(r, Bfunc, μ, t::Number = 0.0; species = Proton)
    return get_adiabaticity(r, Bfunc, species.q, species.m, μ, t)
end

function get_adiabaticity(xv, p::Tuple, t)
    state_len = length(xv)

    if state_len == 6 # FO or Hybrid
        q2m, m, Efunc, Bfunc, _ = p
        q = q2m * m

        x = xv[SA[1, 2, 3]]
        v = xv[SA[4, 5, 6]]

        B = Bfunc(x, t)
        Bmag = norm(B)
        b̂ = B / Bmag
        vpar_val = v ⋅ b̂
        vperp_vec = v - vpar_val * b̂
        μ = m * (vperp_vec ⋅ vperp_vec) / (2 * Bmag)

        return get_adiabaticity(x, Bfunc, q, m, μ, t)
    elseif state_len == 4 # GC
        q, q2m, μ, Efunc, Bfunc = p
        m = q / q2m

        x = xv[SA[1, 2, 3]]
        return get_adiabaticity(x, Bfunc, q, m, μ, t)
    else
        error("Unsupported solution state dimension: $state_len")
    end
end

"""
    get_adiabaticity(sol::AbstractODESolution)

Calculate the adiabaticity parameter `ϵ` from a solution `sol` at all time steps.
The parameters `q`, `m`, `Bfunc` are extracted from `sol.prob.p`.
"""
function get_adiabaticity(sol::AbstractODESolution)
    n = length(sol.t)
    ε = Vector{Float64}(undef, n)
    p = sol.prob.p

    @inbounds for i in 1:n
        ε[i] = get_adiabaticity(sol.u[i], p, sol.t[i])
    end

    return ε
end

"""
    get_adiabaticity(sol::AbstractODESolution, t)

Calculate the adiabaticity parameter `ϵ` from a solution `sol` at time `t`.
`t` must be a scalar.
"""
function get_adiabaticity(sol::AbstractODESolution, t)
    return get_adiabaticity(sol(t), sol.prob.p, t)
end


"""
    get_particle_flux(sol, surface::Union{Disk, Plane, Sphere}; weight=1.0)

Calculate both the velocity flux and the number of particles crossing a virtual detector.
The generic `sol` must store the discrete time steps in `sol.t` and values in `sol.u`.
If `weight` is provided, the result is multiplied by the weight.
Returns a tuple `(velocity_flux, number_flux)`.
For `Disk` and `Sphere`, the returned `velocity_flux` is normalized by the area.
See also [`get_particle_fluxes`](@ref).
"""
function get_particle_flux(sol, surface::Union{Disk, Plane, Sphere}; weight = 1.0)
    t = sol.t
    u = sol.u
    T = eltype(u[1])
    velocities = SVector{3, T}[]
    # Most particles cross a surface once
    sizehint!(velocities, 1)

    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)
    number_flux = zero(weight)

    @inbounds for i in 1:(length(t) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        flux = _check_intersection!(
            velocities, s1, s2, surface, sol, t[i], t[i + 1], weight
        )
        number_flux += flux
        s1 = s2
    end

    _calculate_flux!(velocities, surface)
    return velocities, number_flux
end

"""
    get_particle_fluxes(sols, surface; weights=1.0)::(velocity_fluxes, total_number_flux)

Calculate both the velocity fluxes and total weighted particle flux for an ensemble of trajectories `sols`.
"""
function get_particle_fluxes(sols, surface::Union{Disk, Plane, Sphere}, weights::Number)
    return get_particle_fluxes(sols, surface, Base.Iterators.repeated(weights))
end

function get_particle_fluxes(sols, surface::Union{Disk, Plane, Sphere}, weights)
    T = eltype(first(sols).u[1])
    velocities = SVector{3, T}[]
    sizehint!(velocities, length(sols))
    total_n_flux = zero(eltype(weights))

    @inbounds for (sol, w) in zip(sols, weights)
        t, u = sol.t, sol.u
        u1 = u[1]
        p1 = Point(u1[1], u1[2], u1[3])
        s1 = _signed_distance(p1, surface)

        for i in 1:(length(t) - 1)
            u2 = u[i + 1]
            p2 = Point(u2[1], u2[2], u2[3])
            s2 = _signed_distance(p2, surface)

            flux = _check_intersection!(
                velocities, s1, s2, surface, sol, t[i], t[i + 1], w
            )
            total_n_flux += flux
            s1 = s2
        end
    end

    _calculate_flux!(velocities, surface)

    return velocities, total_n_flux
end

"""
    get_particle_fluxes(sols, surfaces::AbstractVector{<:Union{Disk, Plane, Sphere}}, weights)

Calculate particle fluxes across multiple detectors for an ensemble of trajectories.
"""
function get_particle_fluxes(
        sols, surfaces::AbstractVector{<:Union{Disk, Plane, Sphere}}, weights::Number
    )
    return get_particle_fluxes(sols, surfaces, Base.Iterators.repeated(weights))
end

function get_particle_fluxes(
        sols, surfaces::AbstractVector{<:Union{Disk, Plane, Sphere}}, weights
    )
    nsurfaces = length(surfaces)
    T = eltype(first(sols).u[1])
    results = Vector{Vector{SVector{3, T}}}(undef, nsurfaces)
    @inbounds for j in 1:nsurfaces
        results[j] = sizehint!(SVector{3, T}[], length(sols))
    end
    total_n_fluxes = zeros(eltype(weights), nsurfaces)
    s1_buf = Vector{eltype(weights)}(undef, nsurfaces)

    @inbounds for (sol, w) in zip(sols, weights)
        t, u = sol.t, sol.u
        u1 = u[1]
        p1 = Point(u1[1], u1[2], u1[3])
        for j in eachindex(surfaces)
            s1_buf[j] = _signed_distance(p1, surfaces[j])
        end

        for i in eachindex(t)[1:(end - 1)]
            u2 = u[i + 1]
            p2 = Point(u2[1], u2[2], u2[3])
            for j in eachindex(surfaces)
                s2 = _signed_distance(p2, surfaces[j])
                flux = _check_intersection!(
                    results[j], s1_buf[j], s2, surfaces[j], sol, t[i], t[i + 1], w
                )
                total_n_fluxes[j] += flux
                s1_buf[j] = s2
            end
        end
    end

    @inbounds for j in eachindex(surfaces)
        _calculate_flux!(results[j], surfaces[j])
    end

    return results, total_n_fluxes
end

get_particle_fluxes(sols, surface; weights = 1.0) = get_particle_fluxes(sols, surface, weights)
get_particle_fluxes(sols, surfaces::AbstractVector; weights = 1.0) = get_particle_fluxes(sols, surfaces, weights)

@inline function _check_intersection!(velocities, s1, s2, surface, sol, tl, tr, weight)
    flux = zero(weight)
    if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
        f = s1 / (s1 - s2)
        tc = muladd(f, tr - tl, tl)
        ut = sol(tc)
        xc = Point(ut[1], ut[2], ut[3])

        if _is_valid_intersection(xc, surface)
            vt = SVector(ut[4], ut[5], ut[6])
            push!(velocities, vt * weight)
            flux = weight
        end
    end

    return flux
end

@inline function _signed_distance(p::Point, surface::Disk)
    n = normal(surface.plane)
    v = p - surface.plane.p
    return (v ⋅ n).val
end

@inline function _signed_distance(p::Point, surface::Plane)
    n = normal(surface)
    v = p - surface.p
    return (v ⋅ n).val
end

@inline function _signed_distance(p::Point, surface::Sphere)
    v = p - surface.center
    return (v ⋅ v - surface.radius^2).val
end

@inline function _is_valid_intersection(p::Point, surface::Disk)
    center = surface.plane.p
    return (p - center) ⋅ (p - center) <= surface.radius^2
end

@inline function _is_valid_intersection(p::Point, surface::Union{Plane, Sphere})
    return true
end

function _calculate_flux!(velocities, surface::Union{Disk, Sphere})
    inv_area = inv(area(surface).val)
    @inbounds for i in eachindex(velocities)
        velocities[i] *= inv_area
    end
    return velocities
end

function _calculate_flux!(velocities, ::Plane)
    return velocities
end
