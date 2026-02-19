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
    set_axes_equal(ax)

Set 3D plot axes to equal scale for Matplotlib.
Make axes of 3D plot have equal scale so that spheres appear as spheres and cubes as cubes.
Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don't work on 3D.
"""
function set_axes_equal(ax)
    limits = zeros(2, 3)
    limits[:, 1] .= ax.get_xlim3d()
    limits[:, 2] .= ax.get_ylim3d()
    limits[:, 3] .= ax.get_zlim3d()
    origin = mean(limits, dims = 1)
    radius = @. 0.5 * max(abs(limits[2, :] - limits[1, :]))
    x, y, z = origin
    ax.set_xlim3d([x - radius[1], x + radius[1]])
    ax.set_ylim3d([y - radius[2], y + radius[2]])
    return ax.set_zlim3d([z - radius[3], z + radius[3]])
end

"""
    get_rotation_matrix(axis::AbstractVector, angle) :: SMatrix{3,3}

Create a rotation matrix for rotating a 3D vector around a unit `axis` by an `angle` in
radians.
Reference: [Rotation matrix from axis and angle](https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)

# Example

```julia
using LinearAlgebra
v = [-0.5, 1.0, 1.0]
v̂ = normalize(v)
θ = deg2rad(-74)
R = get_rotation_matrix(v̂, θ)
```
"""
function get_rotation_matrix(v, θ)
    sinθ, cosθ = sincos(θ)
    tmp = 1 - cosθ
    return m = @SMatrix [
        cosθ + v[1]^2 * tmp v[1] * v[2] * tmp - v[3] * sinθ v[1] * v[3] * tmp + v[2] * sinθ;
        v[1] * v[2] * tmp + v[3] * sinθ cosθ + v[2]^2 * tmp v[2] * v[3] * tmp - v[1] * sinθ;
        v[1] * v[3] * tmp - v[2] * sinθ v[3] * v[2] * tmp + v[1] * sinθ cosθ + v[3]^2 * tmp
    ]
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

    v_E = SVector{3}((E × B) / Bmag^2) # ExB drift
    v_gyro = SVector{3}(v - v_E)

    # Calculate perpendicular component
    b̂ = SVector{3}(B / Bmag)
    v_gyro_par = (v_gyro ⋅ b̂) * b̂
    v_gyro_perp = v_gyro - v_gyro_par
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
    m = p[2]
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
    v_E = (E × B) / Bmag^2

    # Determine if relativistic
    func_name = Symbol(sol.prob.f.f)
    is_relativistic = func_name === :trace_relativistic! ||
        func_name === :trace_relativistic ||
        func_name === :trace_relativistic_normalized! ||
        func_name === :trace_relativistic_normalized

    # Calculate physical velocity
    if is_relativistic
        v_real = get_relativistic_v(u)
        v_mag = norm(v_real)
        γ = 1 / sqrt(1 - (v_mag / c)^2)
    else
        v_real = u
        γ = 1.0
    end

    v_gyro = v_real - v_E # Subtract drift velocity

    # Calculate perpendicular component of v_gyro
    b̂ = B / Bmag
    v_gyro_par = (v_gyro ⋅ b̂) * b̂
    v_gyro_perp = v_gyro - v_gyro_par
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
    for is in axes(v, 2)
        γv = @view sol[4:6, is]
        γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
        v² = γ²v² / (1 + γ²v² / c^2)
        γ = 1 / √(1 - v² / c^2)
        for i in axes(v, 1)
            v[i, is] = γv[i] / γ
        end
    end

    return v
end

"""
Return the energy [eV] from relativistic `sol`.
"""
function get_energy(sol::AbstractODESolution; m = mᵢ, q = qᵢ)
    e = Vector{eltype(sol.u[1])}(undef, length(sol))
    for i in eachindex(e)
        γv = @view sol[4:6, i]
        γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
        v² = γ²v² / (1 + γ²v² / c^2)
        γ = 1 / √(1 - v² / c^2)
        e[i] = (γ - 1) * m * c^2 / abs(q)
    end

    return e
end

"""
Calculate the energy [eV] of a relativistic particle from γv in [m/s].
"""
function get_energy(γv; m = mᵢ, q = qᵢ)
    γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
    v² = γ²v² / (1 + γ²v² / c^2)
    γ = 1 / √(1 - v² / c^2)

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

    origin = Tuple(get_val(getproperty(g_min, name)) for name in (:x, :y, :z)[1:dim])

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

        vol = zeros(length(dx), length(dy), length(dz))
        @inbounds for k in eachindex(dz), j in eachindex(dy), i in eachindex(dx)
            vol[i, j, k] = dx[i] * dy[j] * dz[k]
        end
        return vol
    elseif paramdim(grid) == 2
        x, y = xyz
        dx = diff(x)
        dy = diff(y)

        vol = zeros(length(dx), length(dy))
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

    for i in 1:n
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
