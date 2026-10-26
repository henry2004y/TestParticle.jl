# Kept for compatibility until the next breaking release, then removed.

PotentialField(phi, A; kw...) = PotentialField(; phi, A, kw...)
PotentialField(A; phi = nothing, kw...) = PotentialField(; phi, A, kw...)

prepare_canonical(phi, A; kw...) = prepare_canonical(; phi, A, kw...)
prepare_canonical(A; phi = nothing, kw...) = prepare_canonical(; phi, A, kw...)

function get_B_parameters(x, t, Bfunc)
    (; B, Bmag, b̂, ∇B, JB) = get_magnetic_properties(x, t, Bfunc)
    return B, Bmag, b̂, ∇B, JB
end

# The general equations with `q = 1, m = 1` (and `c = 1`) and no external force,
# whatever `p` holds besides `E` and `B`.
_normalized(p) = (; q2m = 1, m = 1, p.E, p.B, F = ZeroField(), c = 1)

# Inlined down to the field evaluation (`get_dv` is `@inline`) so that the unit `q2m`, `m`
# and `c` constant-fold and the rebuilt parameters, which copy the possibly large
# interpolated fields, are never materialized.
trace_normalized(y, p, t) = @inline trace(y, _normalized(p), t)
trace_normalized!(dy, y, p, t) = _store!(dy, trace_normalized(y, p, t))
trace_relativistic_normalized(y, p, t) = @inline trace_relativistic(y, _normalized(p), t)
trace_relativistic_normalized!(dy, y, p, t) =
    _store!(dy, trace_relativistic_normalized(y, p, t))

_lorentz_factor(
    γv,
    ::Union{typeof(trace_relativistic_normalized!), typeof(trace_relativistic_normalized)},
    p
) = sqrt(1 + γv ⋅ γv)
