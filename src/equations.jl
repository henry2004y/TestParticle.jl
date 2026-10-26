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
get_dx!(dx, v, x, p, t) = _store!(dx, v)

"""
    get_dv!(dv, v, x, p, t)

In-place solver components for `DynamicalODEProblem` (velocity).
"""
get_dv!(dv, v, x, p, t) = _store!(dv, @inline get_dv(v, x, p, t))

"""
    get_relativistic_v(γv; c = c)

Velocity of a particle with the relativistic velocity `γv`.
"""
function get_relativistic_v(γv; c = c)
    γ²v² = γv ⋅ γv
    return γv / √(1 + γ²v² / oftype(float(γ²v²), c)^2)
end

# The in-place equations store the out-of-place right-hand side.
@inline function _store!(dy, rhs)
    @inbounds for i in eachindex(rhs)
        dy[i] = rhs[i]
    end
    return
end

"""
    trace!(dy, y, p, t)

ODE equations for charged particle moving in EM field and external force field with in-place form.
"""
trace!(dy, y, p, t) = _store!(dy, @inline trace(y, p, t))

"""
    trace(y, p, t)::SVector{6}

ODE equations for charged particle moving in EM field and external force field with out-of-place form.
"""
function trace(y, p, t)
    v = get_v(y)
    return vcat(v, get_dv(v, y, p, t))
end

_speed_of_light(p) = hasproperty(p, :c) ? p.c : c

"""
    trace_relativistic!(dy, y, p, t)

ODE equations for relativistic charged particle (x, γv) moving in EM field with in-place form.
The speed of light is `p.c` when present, `c` in SI otherwise.
"""
trace_relativistic!(dy, y, p, t) = _store!(dy, @inline trace_relativistic(y, p, t))

"""
    trace_relativistic(y, p, t) -> SVector{6}

ODE equations for relativistic charged particle (x, γv) moving in EM field with out-of-place form.
"""
function trace_relativistic(y, p, t)
    v = get_relativistic_v(get_v(y); c = _speed_of_light(p))
    return vcat(v, get_dv(v, y, p, t))
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

_is_relativistic(f) = f in (trace_relativistic!, trace_relativistic)
