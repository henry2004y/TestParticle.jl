# Guiding center.

"Split a full-orbit state into position and velocity."
@inline _split_xv(xv) = xv[SA[1:3...]], xv[SA[4:6...]]

"""
    _gc_position(x, v, B, q, m, t)

Guiding-center position `X = x - ρ(x)`, with the Larmor radius `ρ` evaluated at
the particle: the zeroth-order map.
"""
@muladd @inline function _gc_position(x, v, B, q, m, t)
    b = B(x, t)
    return x - (b × v) / ((q / m) * (b ⋅ b))
end

@inline _gc_lorentz_factor_full(u, ::Union{GuidingCenter, FullOrbit}, c_val) =
    one(eltype(u))
@inline function _gc_lorentz_factor_full(
        u, ::Union{RelativisticGuidingCenter, RelativisticOrbit}, c_val
    )
    u² = u ⋅ u
    return sqrt(1 + u² / oftype(float(u²), c_val)^2)
end

@inline _gc_to_full_gamma(::GuidingCenter, vpar, w², c_val) = one(w²)
@inline function _gc_to_full_gamma(::RelativisticGuidingCenter, vpar, w², c_val)
    u² = vpar^2 + w²
    return sqrt(1 + u² / oftype(float(u²), c_val)^2)
end

"""
    _gc_frame(x, v, E, B, m, t, model = GuidingCenter(), c_val = c) -> (; b̂, w, vpar, μ)

Field-aligned frame at position `x`: unit field `b̂`, gyration velocity
`w = v - vpar b̂ - γ v_E`, parallel velocity `vpar = b̂ ⋅ v` and magnetic moment
`μ = m|w|² / (2|B|)`. For `RelativisticGuidingCenter`, `v` is the proper velocity `γv_coord`.
"""
@muladd @inline function _gc_frame(
        x, v, E, B, m, t, model = GuidingCenter(), c_val = c
    )
    b = B(x, t)
    Bmag = norm(b)
    inv_Bmag = inv(Bmag)
    b̂ = b * inv_Bmag
    vpar = b̂ ⋅ v
    γ = _gc_lorentz_factor_full(v, model, c_val)
    w = v - vpar * b̂ - γ * ((E(x, t) × b̂) * inv_Bmag)
    return (; b̂, w, vpar, μ = m * (w ⋅ w) * (inv_Bmag / 2))
end

"Gyrophase of the gyration velocity `w` in the frame built from `b̂`."
@inline function _gc_phase(b̂, w)
    e1, e2 = _unit_perp_vectors(b̂)
    return atan(w ⋅ e2, w ⋅ e1)
end

"""
    _get_gc_parameters(xv, E, B, q, m, t, model = GuidingCenter(), c_val = c)

Extract the guiding-center position `X`, parallel (proper) velocity `vpar`, magnetic
moment `μ`, and gyrophase `phase` from a full-particle state `xv` at time `t`.
`E` and `B` must be callable as `E(x, t)` / `B(x, t)`; static fields simply
ignore `t`.

#TODO : support external forces
"""
@inline function _get_gc_parameters(
        xv, E, B, q, m, t, model = GuidingCenter(), c_val = c
    )
    x, v = _split_xv(xv)
    X = _gc_position(x, v, B, q, m, t)
    f = _gc_frame(X, v, E, B, m, t, model, c_val)
    return X, f.vpar, q, q / m, f.μ, _gc_phase(f.b̂, f.w)
end

"""
    _gc_adia_parameters(xv, E, B, q, m, t) -> (X, vpar, μ)

Guiding-center position, parallel velocity, and magnetic moment of the
full-orbit state `xv` — the quantities an adiabaticity estimate needs. The
gyrophase (`get_perp_vector` and an `atan2`) is skipped because the estimate
never uses it; the values match those of [`_get_gc_parameters`](@ref).
"""
@inline function _gc_adia_parameters(xv, E, B, q, m, t)
    x, v = _split_xv(xv)
    X = _gc_position(x, v, B, q, m, t)
    f = _gc_frame(X, v, E, B, m, t)
    return X, f.vpar, f.μ
end

# Guiding-center state `[X..., vpar, μ]` and gyrophase of `xv`, refined by up to `niter`
# fixed-point iterations of `X = x - ρ(X)`: evaluating the Larmor radius at the guiding
# center instead of the particle removes the `O(ε)` error of the zeroth-order map
# (`niter = 0`). The result is the exact inverse of `_gc_to_full`, so a switch to the guiding
# center and back returns the particle where it was. An iteration is kept only if it reduces
# the residual `‖x - (X + ρ(X))‖`, so where the map stops contracting the result is no worse
# than zeroth order.
@inline function _gc_switch_state_checked(xv, E, B, q, m, t, niter)
    x, v = _split_xv(xv)
    X, vpar, _, _, μ, phase = _get_gc_parameters(xv, E, B, q, m, t)
    state = SVector(X[1], X[2], X[3], vpar, μ)
    niter > 0 || return state, phase, true

    x_rec = _gc_to_full(state, E, B, q, m, t, phase)
    residual = norm(x - get_x(x_rec))

    for _ in 1:niter
        X_next = X + (x - get_x(x_rec))
        f = _gc_frame(X_next, v, E, B, m, t)
        phase_next = _gc_phase(f.b̂, f.w)
        state_next = SVector(X_next[1], X_next[2], X_next[3], f.vpar, f.μ)
        x_rec_next = _gc_to_full(state_next, E, B, q, m, t, phase_next)
        residual_next = norm(x - get_x(x_rec_next))
        residual_next < residual || break

        X, state, phase = X_next, state_next, phase_next
        x_rec, residual = x_rec_next, residual_next
    end

    # issue #667: reject the switch when the gyro-orbit spans a field reversal
    # so `X = x - ρ(X)` fails to invert or `B` flips across `X ± ρ`.
    T = typeof(residual)
    ρ_len = norm(x - X)
    tol_pos = 16 * eps(T) * max(norm(x), norm(X), one(T))
    inv_ok = residual <= max(ρ_len / 4, tol_pos)
    B_X = B(X, t)
    B_opp = B(2 * X - x, t)
    B_prod = norm(B_X) * norm(B_opp)
    same_dir = !iszero(B_prod) && (B_X ⋅ B_opp) > (B_prod / 2)
    return state, phase, inv_ok && same_dir
end

@inline function _gc_state_at_switch(xv, E, B, q, m, t, niter)
    state, phase, _ = _gc_switch_state_checked(xv, E, B, q, m, t, niter)
    return state, phase
end

"""
    full_to_gc(xv, p; t = 0, model = GuidingCenter())

Guiding-center state `[X..., vpar, μ]` of the full-orbit state `xv = [x, y, z, vx, vy, vz]`
at time `t` in the fields of `p` from [`prepare`](@ref), in the float type of `xv`; an
`SVector` for a static `xv`. For `model = RelativisticGuidingCenter()`,
`xv = [x..., (γv)...]`, and the state holds `γvpar` and `μ = p⊥² / (2mB)`.

[`TraceGCProblem`](@ref) converts a full-orbit `u0` itself.

# Examples
```jldoctest
julia> B(x) = SA[0.0, 0.0, 1.0e-8];

julia> p = prepare(; B, species = Proton);

julia> u0 = SA[1.0, 0.0, 0.0, 0.0, 1.0, 0.0];

julia> sol = solve(TraceGCProblem(u0, (0.0, 1.0), p), Vern7());

julia> isapprox(sol.u[end], TestParticle.full_to_gc(u0, p); atol = 1.0e-12)
true
```
"""
function full_to_gc(xv, p; t = 0, model = GuidingCenter())
    T = float(eltype(xv))
    X, vpar, _, _, μ, _ = _get_gc_parameters(
        xv, p.E, p.B, p.q, p.m, t, model, _speed_of_light(p)
    )
    state = SVector{5, T}(X[1], X[2], X[3], vpar, μ)
    return xv isa StaticArray ? state : Vector(state)
end

"""
    gc_to_full(state_gc, p; phase = 0, t = 0, model = GuidingCenter())

Full-orbit state `[x, y, z, vx, vy, vz]` (`[x..., (γv)...]` for
`RelativisticGuidingCenter()`) of the guiding-center state `[X..., vpar, μ]` at gyrophase
`phase` and time `t`: the inverse of [`full_to_gc`](@ref).
"""
function gc_to_full(state_gc, p; phase = 0, t = 0, model = GuidingCenter())
    length(state_gc) == 5 || throw(
        ArgumentError(
            "A guiding-center state is `[X..., v∥, μ]`, not of length $(length(state_gc))."
        )
    )
    return Vector(
        _gc_to_full(
            state_gc, p.E, p.B, p.q, p.m, t, phase, model,
            _speed_of_light(p)
        )
    )
end

# Takes the fields and charges instead of `p` so the hybrid solver can share it, and returns
# an `SVector{6}` to stay allocation-free in its loop; `gc_to_full` wraps it.
@inline function _gc_to_full(
        state_gc, E_field, B_field, q, m, t, phase = 0,
        model = GuidingCenter(), c_val = c
    )
    R = get_x(state_gc)
    B = B_field(R, t)
    return _gc_to_full(state_gc, E_field, B, norm(B), q, m, t, phase, model, c_val)
end

@muladd @inline function _gc_to_full(
        state_gc, E_field, B::SVector{3}, Bmag::Real, q, m, t, phase = 0,
        model = GuidingCenter(), c_val = c
    )
    R = get_x(state_gc)
    vpar, μ = state_gc[4], state_gc[5]

    E = E_field(R, t)
    inv_Bmag = inv(Bmag)
    b̂ = B * inv_Bmag

    v_E = (E × b̂) * inv_Bmag

    # perp speed, μ = m * w^2 / (2B)
    w² = 2 * (μ / m) * Bmag
    w = sqrt(w²)
    γ = _gc_to_full_gamma(model, vpar, w², c_val)

    # perp vector
    e1, e2 = _unit_perp_vectors(b̂)
    s, c_ph = sincos(phase)
    v_gyr = w * (c_ph * e1 + s * e2)
    v_perp = v_gyr + γ * v_E

    v = vpar * b̂ + v_perp

    # gyroradius
    inv_Ω = (m / q) * inv_Bmag
    ρ_vec = (b̂ × v_perp) * inv_Ω

    x = R + ρ_vec

    return SVector{6}(x[1], x[2], x[3], v[1], v[2], v[3])
end

"""
    get_gc(xu, param; use_vE = false, model = GuidingCenter())
    get_gc(x, y, z, vx, vy, vz, bx, by, bz, q2m)

Calculate the coordinates of the guiding center according to the phase space coordinates of a particle.
Reference: [wiki](https://en.wikipedia.org/wiki/Guiding_center)

Nonrelativistic definition:

```math
\\mathbf{X}=\\mathbf{x}-m\\frac{\\mathbf{b}\\times\\mathbf{v}}{qB}
```

For relativistic states `xu = [x..., (γv)...]`, `ρ = (b × γv) / (qB / m)` has the same form
when `use_vE = false`; with `use_vE = true` and `model = RelativisticGuidingCenter()` (or
`RelativisticOrbit()`), the drift subtraction uses `γv - γv_E`.

If `use_vE=true`, it calculates the guiding center relative to the E × B drift frame:

```math
\\mathbf{X}=\\mathbf{x}-m\\frac{\\mathbf{b}\\times(\\mathbf{v}-\\mathbf{v}_E)}{qB}
```

which is essential for correctly identifying the polarization drift in time-varying fields.
"""
function get_gc(xu, param; use_vE::Bool = false, model = GuidingCenter())
    # TODO: support external forces
    q2m = param.q2m
    B_field = param.B
    t = length(xu) == 7 ? xu[end] : zero(eltype(xu))
    v = xu[SA[4:6...]]
    B = B_field(xu, t)
    B² = B[1]^2 + B[2]^2 + B[3]^2

    if use_vE
        E_field = param.E
        E = E_field(xu, t)
        γ = _gc_lorentz_factor_full(v, model, _speed_of_light(param))
        v = v - γ * ((E × B) / B²)
    end

    # vector of Larmor radius
    ρ = B × v ./ (q2m * B²)

    return X = @views xu[1:3] - ρ
end

function get_gc(x, y, z, vx, vy, vz, bx, by, bz, q2m)
    l = SVector{3}(x, y, z)
    B = SVector{3}(bx, by, bz)
    v = SVector{3}(vx, vy, vz)

    B² = bx^2 + by^2 + bz^2
    # vector of Larmor radius
    ρ = B × v ./ (q2m * B²)

    return X = l - ρ
end
