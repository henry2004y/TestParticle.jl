# Tracing equations.

get_x(u) = @inbounds SA[u[1], u[2], u[3]]
get_v(u) = @inbounds SA[u[4], u[5], u[6]]

@inline function get_dv(v::AbstractVector{T}, x, p, t) where {T <: AbstractFloat}
    E = _eval_field(T, p.E, x, t)
    B = _eval_field(T, p.B, x, t)
    F = _eval_field(T, p.F, x, t)

    return SVector{3, T}(T(p.q2m) * (v × B + E) + F / T(p.m))
end

@inline get_dv(v, x, p, t) = p.q2m * (v × p.B(x, t) + p.E(x, t)) + p.F(x, t) / p.m

"""
    get_dx!(dx, v, x, p, t)

In-place solver components for `DynamicalODEProblem` (location).
"""
function get_dx!(dx, v, x, p, t)
    @inbounds for i in eachindex(dx, v)
        dx[i] = v[i]
    end
    return
end

"""
    get_dv!(dv, v, x, p, t)

In-place solver components for `DynamicalODEProblem` (velocity).
"""
function get_dv!(dv, v, x, p, t)
    val = get_dv(v, x, p, t)
    @inbounds for i in eachindex(dv, val)
        dv[i] = val[i]
    end

    return
end

function get_relativistic_v(γv::AbstractVector{T}; c = c) where {T <: AbstractFloat}
    c_val = T(c)
    γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
    if γ²v² > eps(T)
        v̂ = normalize(γv)
    else # no velocity
        v̂ = SVector{3, T}(0, 0, 0)
    end
    return √(γ²v² / (1 + γ²v² / c_val^2)) * v̂
end

function get_relativistic_v(γv; c = c)
    γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
    if γ²v² > eps(eltype(γv))
        v̂ = normalize(γv)
    else # no velocity
        v̂ = zero(γv)
    end
    return √(γ²v² / (1 + γ²v² / c^2)) * v̂
end

"""
    trace!(dy, y, p, t)

ODE equations for charged particle moving in EM field and external force field with in-place form.
"""
function trace!(dy, y, p, t)
    v = get_v(y)
    dv = get_dv(v, y, p, t)
    @inbounds dy[1] = v[1]
    @inbounds dy[2] = v[2]
    @inbounds dy[3] = v[3]
    @inbounds dy[4] = dv[1]
    @inbounds dy[5] = dv[2]
    @inbounds dy[6] = dv[3]

    return
end

"""
    trace(y, p, t)::SVector{6}

ODE equations for charged particle moving in EM field and external force field with out-of-place form.
"""
function trace(y, p, t)
    v = y[SA[4:6...]]
    dv = get_dv(v, y, p, t)
    return vcat(v, dv)
end

_speed_of_light(p) = hasproperty(p, :c) ? p.c : c

"""
    trace_relativistic!(dy, y, p, t)

ODE equations for relativistic charged particle (x, γv) moving in EM field with in-place form.
The speed of light is `p.c` when present, `c` in SI otherwise.
"""
function trace_relativistic!(dy, y, p, t)
    γv = get_v(y)
    v = get_relativistic_v(γv; c = _speed_of_light(p))
    @inbounds dy[1:3] = v
    @inbounds dy[4:6] = get_dv(v, y, p, t)

    return
end

"""
    trace_relativistic(y, p, t) -> SVector{6}

ODE equations for relativistic charged particle (x, γv) moving in EM field with out-of-place form.
"""
function trace_relativistic(y, p, t)
    γv = get_v(y)
    v = get_relativistic_v(γv; c = _speed_of_light(p))
    dv = get_dv(v, y, p, t)

    return vcat(v, dv)
end

"""
    _lorentz_factor(u, f, p)

Lorentz factor of the velocity part `u` of a state traced with the equations `f` and
parameters `p`; the relativistic equations store `γv` in place of `v`.
"""
_lorentz_factor(u, f, p) = one(eltype(u))
_lorentz_factor(
    γv, ::Union{typeof(trace_relativistic!), typeof(trace_relativistic)}, p
) = sqrt(1 + (γv ⋅ γv) / _speed_of_light(p)^2)

# Kept for compatibility until the next breaking release: the general equations with
# `q = 1, m = 1` (and `c = 1`) and no external force, whatever `p` holds besides `E` and `B`.
_normalized(p) = (; q2m = 1, m = 1, p.E, p.B, F = ZeroField(), c = 1)

# Inlined down to the field evaluation (`get_dv` is `@inline`) so that the unit `q2m`, `m`
# and `c` constant-fold and the rebuilt parameters, which copy the possibly large
# interpolated fields, are never materialized.
trace_normalized!(dy, y, p, t) = @inline trace!(dy, y, _normalized(p), t)
trace_normalized(y, p, t) = @inline trace(y, _normalized(p), t)
trace_relativistic_normalized!(dy, y, p, t) =
    @inline trace_relativistic!(dy, y, _normalized(p), t)
trace_relativistic_normalized(y, p, t) = @inline trace_relativistic(y, _normalized(p), t)

_lorentz_factor(
    γv,
    ::Union{typeof(trace_relativistic_normalized!), typeof(trace_relativistic_normalized)},
    p
) = sqrt(1 + γv ⋅ γv)
