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

"""
    get_number_density_flux(grid::CartesianGrid, sols, dt)

Calculate the steady state particle number density flux on a uniform Cartesian grid.
The flux is estimated by accumulating the number of particles in each cell at time steps `dt`,
divided by the surface area of each cell.

# Arguments

  - `grid`: A `CartesianGrid` from `Meshes.jl`.
  - `sols`: Particle trajectory solutions (e.g. `EnsembleSolution`).
  - `dt`: Time step for sampling particle positions.
"""
function get_number_density_flux(grid::CartesianGrid, sols, dt)
    counts = zeros(Int, size(grid))
    dim = paramdim(grid)

    g_min = coords(minimum(grid))
    Δx = spacing(grid)

    get_val(x) = hasproperty(x, :val) ? x.val : x

    origin = if dim == 3
        (get_val(g_min.x), get_val(g_min.y), get_val(g_min.z))
    elseif dim == 2
        (get_val(g_min.x), get_val(g_min.y))
    else
        (get_val(g_min.x),)
    end

    spacings = Tuple(get_val(d) for d in Δx)

    if dim == 3
        dx, dy, dz = spacings
        cell_area = 2 * (dx * dy + dy * dz + dz * dx)
    elseif dim == 2
        dx, dy = spacings
        cell_area = dx * dy
    elseif dim == 1
        cell_area = 1.0
    end

    sz = size(grid)

    for sol in sols
        t0 = sol.prob.tspan[1]
        t1 = sol.prob.tspan[2]
        t_start, t_end = t0 < t1 ? (t0, t1) : (t1, t0)

        for t in t_start:dt:t_end
            state = sol(t)
            pos = state[1:dim]

            idx = MVector{dim, Int}(undef)
            in_bounds = true

            for d in 1:dim
                val = floor(Int, (pos[d] - origin[d]) / spacings[d])
                id = val + 1
                if 1 <= id <= sz[d]
                    idx[d] = id
                else
                    in_bounds = false
                    break
                end
            end

            if in_bounds
                counts[CartesianIndex(Tuple(idx))] += 1
            end
        end
    end

    return counts ./ cell_area
end

"""
    get_number_density(sols, grid, t)

Calculate particle number density at time `t` in a given `grid`.
"""
function get_number_density(sols, grid, t)
    dims = size(grid)
    counts = zeros(Int, dims)
    ranges = _get_ranges_val(makegrid(grid))
    dim = paramdim(grid)

    for sol in sols
        # Check if the solution is defined at time t
        if sol.t[1] <= t <= sol.t[end]
            pos = _strip_units(sol(t))

            indices = ntuple(i -> searchsortedlast(ranges[i], pos[i]), dim)
            inbounds = all(ntuple(i -> 1 <= indices[i] <= dims[i], dim))
            if inbounds
                counts[indices...] += 1
            end
        end
    end

    # Divide by cell volume to get density
    vol = _get_cell_volume(grid)

    return counts ./ vol
end

"""
    get_number_density(sols, grid, t_start, t_end, dt)

Calculate time-averaged particle number density from `t_start` to `t_end` with step `dt`.
"""
function get_number_density(sols, grid, t_start, t_end, dt)
    dims = size(grid)
    counts = zeros(Int, dims)
    ranges = _get_ranges_val(makegrid(grid))
    dim = paramdim(grid)

    # For averaging, we accumulate counts at each time step
    t_steps = t_start:dt:t_end
    n_steps = length(t_steps)

    for t in t_steps
        for sol in sols
            if sol.t[1] <= t <= sol.t[end]
                pos = _strip_units(sol(t))

                indices = ntuple(i -> searchsortedlast(ranges[i], pos[i]), dim)
                inbounds = all(ntuple(i -> 1 <= indices[i] <= dims[i], dim))
                if inbounds
                    counts[indices...] += 1
                end
            end
        end
    end

    vol = _get_cell_volume(grid)

    return (counts ./ n_steps) ./ vol
end

"Helper function to strip units from values."
_strip_units(x::AbstractArray) = map(v -> _strip_units(v), x)

_strip_units(x) =
if hasproperty(x, :val)
    return x.val
else
    return x
end

_get_ranges_val(ranges::Tuple) = map(_strip_units, ranges)

function _get_cell_volume(grid::CartesianGrid)
    sp = spacing(grid)
    # Strip units from spacing if necessary, though spacing typically returns Quantities if grid is unitful.
    # If we want pure float volume, we need to strip.
    dx = _strip_units(sp[1])
    dy = _strip_units(sp[2])

    if paramdim(grid) == 3
        dz = _strip_units(sp[3])
        return dx * dy * dz
    elseif paramdim(grid) == 2
        return dx * dy
    end
end

function _get_cell_volume(grid::RectilinearGrid)
    # Variable cell volume
    # This is tricky because we return a single array of densities, but cells have different volumes.
    # We should return density per cell, so we divide counts[i,j,k] by volume[i,j,k].

    xyz = makegrid(grid) # Returns ranges or vectors

    # Strip units for volume calculation if needed, or keep them if consistent.
    xyz = _get_ranges_val(xyz)

    if paramdim(grid) == 3
        x, y, z = xyz
        dx = diff(x)
        dy = diff(y)
        dz = diff(z)

        vol = Array{eltype(dx), 3}(undef, length(dx), length(dy), length(dz))
        @inbounds for k in eachindex(dz), j in eachindex(dy), i in eachindex(dx)
            vol[i, j, k] = dx[i] * dy[j] * dz[k]
        end
        return vol
    elseif paramdim(grid) == 2
        x, y = xyz
        dx = diff(x)
        dy = diff(y)

        vol = Array{eltype(dx), 2}(undef, length(dx), length(dy))
        @inbounds for j in eachindex(dy), i in eachindex(dx)
            vol[i, j] = dx[i] * dy[j]
        end
        return vol
    end
end

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
    T = eltype(sol.u[1])
    velocities = SVector{3, T}[]

    u1 = sol.u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)

    @inbounds for i in eachindex(t)[1:(end - 1)]
        u2 = sol.u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        # Check if the line segment intersects the surface
        if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
            f = s1 / (s1 - s2)
            tc = t[i] + f * (t[i + 1] - t[i])
            ut = sol(tc)
            xc = Point(ut[1], ut[2], ut[3])

            if _is_valid_intersection(xc, surface)
                push!(velocities, SVector{3, T}(ut[4], ut[5], ut[6]))
            end
        end
        u1, p1, s1 = u2, p2, s2
    end

    return _calculate_flux(velocities, surface, weight)
end

"""
    get_particle_fluxes(sols, surface; weights=1.0)::(velocity_fluxes, total_number_flux)

Calculate both the velocity fluxes and total weighted particle flux for an ensemble of trajectories `sols`.
"""
function get_particle_fluxes(sols, surface, weights::Number)
    T = eltype(first(sols).u[1])
    velocities = SVector{3, T}[]
    total_n_flux = 0.0
    for sol in sols
        v_f, n_f = get_particle_flux(sol, surface; weight = weights)
        append!(velocities, v_f)
        total_n_flux += n_f
    end
    return velocities, total_n_flux
end

function get_particle_fluxes(sols, surface, weights::AbstractVector)
    T = eltype(first(sols).u[1])
    velocities = SVector{3, T}[]
    total_n_flux = 0.0
    for (sol, w) in zip(sols, weights)
        v_f, n_f = get_particle_flux(sol, surface; weight = w)
        append!(velocities, v_f)
        total_n_flux += n_f
    end
    return velocities, total_n_flux
end

get_particle_fluxes(sols, surface; weights = 1.0) = get_particle_fluxes(sols, surface, weights)

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

function _calculate_flux(velocities, surface::Disk, weight)
    area = π * surface.radius^2
    n_flux = length(velocities) * weight
    v_flux = (velocities .* weight) ./ area.val
    return v_flux, n_flux
end

function _calculate_flux(velocities, surface::Plane, weight)
    n_flux = length(velocities) * weight
    v_flux = velocities .* weight
    return v_flux, n_flux
end

function _calculate_flux(velocities, surface::Sphere, weight)
    area = 4π * surface.radius^2
    n_flux = length(velocities) * weight
    v_flux = (velocities .* weight) ./ area.val
    return v_flux, n_flux
end
