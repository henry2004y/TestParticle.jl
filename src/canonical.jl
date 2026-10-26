# Canonical equations.

# ∇(A ⋅ v) at fixed v; dropped by dispatch for an absent `A`, as `_eval_field` does.
@inline _grad_A_v(grad_A, x, t, v) = grad_A(x, t) * v
_grad_A_v(::ZeroField, x, t, v) = ZeroVector()

@inline _kinetic_momentum(x, p_can, p, t) =
    p_can - p.q * _eval_field(eltype(x), p.pf.A, x, t)

@inline _kinetic_to_velocity(p_kin, p, relativistic::Bool) = relativistic ?
    p_kin / (p.m * √(1 + (p_kin ⋅ p_kin) / (p.m * p.c)^2)) : p_kin / p.m

@inline function _trace_canonical(y, p, t, relativistic::Bool)
    (; q, pf) = p
    x = get_x(y)
    p_can = get_v(y)
    v = _kinetic_to_velocity(_kinetic_momentum(x, p_can, p, t), p, relativistic)

    grad_phi = _eval_field(eltype(x), pf.grad_phi, x, t)
    dp = SVector{3}(q * _grad_A_v(pf.grad_A, x, t, v) - q * grad_phi)
    return vcat(v, dp)
end

trace_canonical(y, p, t) = _trace_canonical(y, p, t, false)
trace_canonical!(dy, y, p, t) = _store!(dy, @inline trace_canonical(y, p, t))
trace_canonical_relativistic(y, p, t) = _trace_canonical(y, p, t, true)
trace_canonical_relativistic!(dy, y, p, t) =
    _store!(dy, @inline trace_canonical_relativistic(y, p, t))

function velocity_to_canonical(
        x0::AbstractVector, v0::AbstractVector, p;
        relativistic::Bool = false, t = 0.0
    )
    (; q, m, c, pf) = p
    x = SA[x0[1], x0[2], x0[3]]
    v = SA[v0[1], v0[2], v0[3]]
    γ = relativistic ? 1 / √(1 - (v ⋅ v) / c^2) : 1
    return vcat(x, γ * m * v + q * _eval_field(eltype(x), pf.A, x, t))
end

function canonical_to_velocity(
        x::AbstractVector, p_can::AbstractVector, p;
        relativistic::Bool = false, t = 0.0
    )
    x_vec = SA[x[1], x[2], x[3]]
    p_vec = SA[p_can[1], p_can[2], p_can[3]]
    return _kinetic_to_velocity(_kinetic_momentum(x_vec, p_vec, p, t), p, relativistic)
end

function canonical_to_velocity(u::AbstractVector, p; relativistic::Bool = false, t = 0.0)
    return canonical_to_velocity(get_x(u), get_v(u), p; relativistic, t)
end

function canonical_hamiltonian(u::AbstractVector, p, t = 0.0; relativistic::Bool = false)
    (; q, m, c, pf) = p
    x = get_x(u)
    p_kin = _kinetic_momentum(x, get_v(u), p, t)
    p² = p_kin ⋅ p_kin
    K = relativistic ? √(m^2 * c^4 + c^2 * p²) : p² / (2 * m)
    phi = pf.phi isa ZeroField ? zero(eltype(x)) : pf.phi(x, t)
    return K + q * phi
end
