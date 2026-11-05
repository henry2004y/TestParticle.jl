# Collection of utility functions and commonly used constants.

include("constants.jl")
include("zero.jl")
include("phasespace.jl")

# Meshes.jl grid and virtual detector function stubs implemented in package extensions

"""
    makegrid(grid::Union{CartesianGrid, RectilinearGrid})

Return the coordinate ranges of a Meshes.jl grid, one range per dimension, e.g.
`(x, y, z)` for a 3-D grid. Requires Meshes.jl to be loaded.
"""
function makegrid end

"""
    get_cell_centers(grid::Union{CartesianGrid, RectilinearGrid})

Return the cell center coordinates of a Meshes.jl grid, one range per dimension.
Requires Meshes.jl to be loaded.
"""
function get_cell_centers end

"""
    get_particle_crossings(sol, surface, weight = 1.0) -> (velocities, weights)
    get_particle_crossings(sols, surface, weights = 1.0) -> (velocities, weights)
    get_particle_crossings(sols, surfaces, weights = 1.0) -> (velocities, weights)

Detect every crossing of a virtual detector `surface` and return the velocity [m/s] at
each crossing together with the corresponding weight.

A trajectory can cross the same detector more than once, so a crossing carries the weight
of its trajectory rather than of the trajectory set: passing a scalar `weight` applies it
to every trajectory of the ensemble. With a vector of `surfaces` the two returned vectors
hold one vector per surface. Requires Meshes.jl to be loaded.
"""
function get_particle_crossings end

"""
    get_first_crossing(sol, surface) -> SVector{6}
    get_first_crossing(sols, surface) -> Vector{SVector{6}}

State `[x, y, z, vx, vy, vz]` at the *first* crossing of `surface`, linearly interpolated
between the two bracketing saved points. Returns a state filled with `NaN` for trajectories
that never reach `surface`, which is how unreachable velocities are flagged when tracing a
detector grid backward to a source. Requires Meshes.jl to be loaded.
"""
function get_first_crossing end

"""
    get_particle_flux(sol, surface, weight = 1.0) -> (number_flux, velocity_flux)
    get_particle_flux(sols, surface, weights = 1.0) -> (number_flux, velocity_flux)
    get_particle_flux(sols, surfaces, weights = 1.0) -> (number_fluxes, velocity_fluxes)

Number flux density and velocity flux density through `surface`, or through each of
`surfaces`, obtained by summing the crossing weights of a single trajectory or of an
ensemble of trajectories `sols` and dividing by the surface area. Requires Meshes.jl to be
loaded.
"""
function get_particle_flux end

"""
    _rng_seed(seed::Integer, i::Integer)

Combine the ensemble `seed` and the trajectory index `i` into a well-mixed `UInt64`
used to seed a per-trajectory RNG.

A naive `seed + i` couples the two linearly: incrementing `seed` merely *rotates*
which RNG each trajectory index draws from, so two nearby seeds (e.g. 1 and 2)
produce ensembles that share almost all trajectories and are therefore not
independent. The SplitMix64-style finalizer below is a bijection on `UInt64`, so
any distinct `seed` decorrelates *every* trajectory while remaining deterministic
and reproducible (no reliance on `hash`, which is salted per Julia session).

Note: Julia's multi-argument `Xoshiro(seed, i, …)` constructor does *not* mix —
it embeds the integers directly into the 4×`UInt64` state (`Xoshiro(s0,s1,s2,s3)`
→ `initstate!`). Because `xoshiro256++`'s first output depends only on `s0` and
`s3`, seeding with `Xoshiro(seed, i, 0, 0)` yields the *same* first draw for every
trajectory index `i` and does not decorrelate adjacent seeds. Robust mixing in
`Xoshiro` only happens via the single-argument constructor (SHA2-256 `hash_seed`),
so `seed` and `i` must be folded into one well-mixed integer first — which is what
this function does.
"""
function _rng_seed(seed::Integer, i::Integer)
    z = (seed % UInt64) + 0x9e3779b97f4a7c15 * (i % UInt64)
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    z = z ⊻ (z >> 31)
    return z
end

"""
    _foreach_chunk(f, ensemblealg, n)

Call `f(irange)` on chunks of the trajectory indices `1:n`, all at once for
`EnsembleSerial` and one chunk per thread for `EnsembleThreads`.
"""
_foreach_chunk(f, ::EnsembleSerial, n) = f(1:n)

function _foreach_chunk(f, ::EnsembleThreads, n)
    Threads.@threads for irange in index_chunks(1:n; n = Threads.nthreads())
        f(irange)
    end
    return
end

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
    ϕ = atan(y, x)
    # Not `mod2pi`: it computes in Float64, which device kernels may lack or throttle.
    ϕ < 0 && (ϕ += 2 * oftype(ϕ, π))
    return SVector{3}(r, θ, ϕ)
end

@inline @inbounds cart2sph(x) = cart2sph(x[1], x[2], x[3])

"""
    sph2cartvec(vr, vθ, vϕ, θ, ϕ)

Convert a vector from spherical to Cartesian coordinates. `θ` and `ϕ` define the position
(polar and azimuthal angles in radians) where the spherical components `(vr, vθ, vϕ)` are
defined.
"""
function sph2cartvec(vr, vθ, vϕ, θ, ϕ)
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    vx = sinθ * cosϕ * vr + cosθ * cosϕ * vθ - sinϕ * vϕ
    vy = sinθ * sinϕ * vr + cosθ * sinϕ * vθ + cosϕ * vϕ
    vz = cosθ * vr - sinθ * vθ
    return SVector{3}(vx, vy, vz)
end

"""
    cart2sphvec(vx, vy, vz, x, y, z)

Convert a vector from Cartesian to spherical coordinates at the Cartesian position `(x, y, z)`.
"""
function cart2sphvec(vx, vy, vz, x, y, z)
    ρ = hypot(x, y)
    r = hypot(ρ, z)
    if r == 0
        return SVector{3, eltype(r)}(0, 0, 0)
    end
    vr = (x * vx + y * vy + z * vz) / r
    if ρ == 0
        vθ = z > 0 ? vx : -vx
        vϕ = vy
    else
        vθ = (z * (x * vx + y * vy) / ρ - ρ * vz) / r
        vϕ = (x * vy - y * vx) / ρ
    end
    return SVector{3}(vr, vθ, vϕ)
end

@inline function _unit_perp_vectors(b̂::SVector{3, T}) where {T}
    z = zero(T)
    e1 = if abs(b̂[3]) < T(0.9)
        inv_n = inv(sqrt(muladd(b̂[1], b̂[1], b̂[2] * b̂[2])))
        SA[-b̂[2] * inv_n, b̂[1] * inv_n, z]
    else
        inv_n = inv(sqrt(muladd(b̂[1], b̂[1], b̂[3] * b̂[3])))
        SA[b̂[3] * inv_n, z, -b̂[1] * inv_n]
    end
    return e1, b̂ × e1
end

"""
    get_perp_vector(b::AbstractVector)

Obtain two unit vectors `e1` and `e2` such that `(e1, e2, b)` form a right-handed orthonormal
system.
"""
function get_perp_vector(b)
    b̂ = normalize(SVector{3}(b[1], b[2], b[3]))
    return _unit_perp_vectors(b̂)
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
    model = _model(sol.prob)
    # Interpolate state at time t
    xu = sol(t)
    x = get_x(xu)
    u = get_v(xu)

    p = sol.prob.p
    q2m = p.q2m
    E = p.E(x, t)
    B = p.B(x, t)
    Bmag = norm(B)
    if Bmag == 0
        return Inf
    end

    # Calculate ExB drift
    v_E = _get_exb_drift(E, B, Bmag)

    γ = _lorentz_factor(u, model, p)

    v_gyro = u / γ - v_E # Subtract drift velocity

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
    get_velocity(sol) -> Matrix

Velocity of a full-orbit `sol` at every saved state, as a `3 × n` matrix.
"""
function get_velocity(sol)
    model = _model(sol.prob)
    p = sol.prob.p
    v = Array{eltype(sol.u[1]), 2}(undef, 3, length(sol.u))
    @inbounds for is in axes(v, 2)
        u = get_v(sol.u[is])
        v[:, is] = u / _lorentz_factor(u, model, p)
    end

    return v
end

"""
    get_energy(sol::AbstractODESolution) -> Vector

Kinetic energy per unit charge [eV for `|q| = e`] at every saved state of `sol`, from the
model and the particle of `sol.prob`.
"""
function get_energy(sol::AbstractODESolution)
    p = sol.prob.p
    model = _model(sol.prob)
    return map((u, t) -> _energy(model, u, p, t), sol.u, sol.t)
end

_energy(::FullOrbit, u, p, t) = (v = get_v(u); (v ⋅ v) / (2 * abs(p.q2m)))

function _energy(model::RelativisticOrbit, u, p, t)
    γv = get_v(u)
    return (γv ⋅ γv) / ((_lorentz_factor(γv, model, p) + 1) * abs(p.q2m))
end

function _energy(::GuidingCenter, u, p, t)
    return u[4]^2 / (2 * abs(p.q2m)) + p.μ * norm(p.B(get_x(u), t)) / abs(p.q)
end

"""
Return velocity magnitude from energy in [eV].
"""
energy2velocity(Ek; m = mᵢ, q = qᵢ) = c * sqrt(1 - 1 / (1 + Ek * abs(q) / (m * c^2))^2)

"""
    sample_unit_sphere(rng::AbstractRNG=default_rng())

Sample a unit vector on a sphere uniformly.
"""
function sample_unit_sphere(rng::AbstractRNG = default_rng())
    ϕ = 2π * rand(rng)
    cosθ = 2 * rand(rng) - 1
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
    sample_maxwellian([rng], Tn, m; offset=0.0, u0=SA[0.0, 0.0, 0.0])
    sample_maxwellian([rng], Tn, species::String; offset=0.0, u0=SA[0.0, 0.0, 0.0])

Sample a velocity [m/s] from a Maxwellian distribution with temperature `Tn` [K], mass `m`
[kg], energy offset `offset` [J], and bulk velocity `u0` [m/s].
This function requires `VelocityDistributionFunctions.jl` to be loaded.
"""
function sample_maxwellian end

"""
    TerminateOutside(isoutofdomain)

Create a `DiscreteCallback` that terminates the simulation if
`isoutofdomain(u, p, t)` is true. The condition follows the same `(u, p, t)`
convention as the `isoutside` keyword, and is adapted here to the
`(u, t, integrator)` signature that SciML callbacks use.
"""
function TerminateOutside(isoutofdomain)
    condition = (u, t, integrator) -> isoutofdomain(u, integrator.p, t)
    return DiscreteCallback(condition, terminate!)
end

"""
    get_fields(sol::AbstractODESolution)

Return the electric and magnetic fields at the saved states of `sol`, in the element type
of the state.
"""
function get_fields(sol::AbstractODESolution)
    p = sol.prob.p
    Efunc, Bfunc = p.E, p.B

    E = map((u, t) -> SVector{3, eltype(u)}(Efunc(get_x(u), t)), sol.u, sol.t)
    B = map((u, t) -> SVector{3, eltype(u)}(Bfunc(get_x(u), t)), sol.u, sol.t)

    return E, B
end
