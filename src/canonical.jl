# Canonical equations.

# ∇(A ⋅ v) at fixed v; dropped by dispatch for an absent `A`, as `_eval_field` does.
@inline _grad_A_v(grad_A, x, t, v) = grad_A(x, t) * v
_grad_A_v(::ZeroField, x, t, v) = ZeroVector()

function trace_canonical!(dy, y, p, t)
    (; q, m, pf) = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = _eval_field(eltype(x), pf.A, x, t)
    v = (p_can - q * A) / m

    grad_phi = _eval_field(eltype(x), pf.grad_phi, x, t)
    dp = SVector{3}(q * _grad_A_v(pf.grad_A, x, t, v) - q * grad_phi)

    @inbounds dy[1] = v[1]
    @inbounds dy[2] = v[2]
    @inbounds dy[3] = v[3]
    @inbounds dy[4] = dp[1]
    @inbounds dy[5] = dp[2]
    @inbounds dy[6] = dp[3]
    return
end

function trace_canonical(y, p, t)
    (; q, m, pf) = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = _eval_field(eltype(x), pf.A, x, t)
    v = (p_can - q * A) / m

    grad_phi = _eval_field(eltype(x), pf.grad_phi, x, t)
    dp = SVector{3}(q * _grad_A_v(pf.grad_A, x, t, v) - q * grad_phi)
    return vcat(v, dp)
end

function trace_canonical_relativistic!(dy, y, p, t)
    (; q, m, c, pf) = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = _eval_field(eltype(x), pf.A, x, t)
    p_kin = p_can - q * A
    γ = √(1 + sum(p_kin .^ 2) / (m^2 * c^2))
    v = p_kin / (γ * m)

    grad_phi = _eval_field(eltype(x), pf.grad_phi, x, t)
    dp = SVector{3}(q * _grad_A_v(pf.grad_A, x, t, v) - q * grad_phi)

    @inbounds dy[1] = v[1]
    @inbounds dy[2] = v[2]
    @inbounds dy[3] = v[3]
    @inbounds dy[4] = dp[1]
    @inbounds dy[5] = dp[2]
    @inbounds dy[6] = dp[3]
    return
end

function trace_canonical_relativistic(y, p, t)
    (; q, m, c, pf) = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = _eval_field(eltype(x), pf.A, x, t)
    p_kin = p_can - q * A
    γ = √(1 + sum(p_kin .^ 2) / (m^2 * c^2))
    v = p_kin / (γ * m)

    grad_phi = _eval_field(eltype(x), pf.grad_phi, x, t)
    dp = SVector{3}(q * _grad_A_v(pf.grad_A, x, t, v) - q * grad_phi)
    return vcat(v, dp)
end

function velocity_to_canonical(
        x0::AbstractVector, v0::AbstractVector, p;
        relativistic::Bool = false, t = 0.0
    )
    (; q, m, c, pf) = p
    x_vec = SA[x0[1], x0[2], x0[3]]
    v_vec = SA[v0[1], v0[2], v0[3]]
    A = _eval_field(eltype(x_vec), pf.A, x_vec, t)

    p_can = if relativistic
        v2 = sum(v_vec .^ 2)
        γ = 1 / √(1 - v2 / c^2)
        γ * m * v_vec + q * A
    else
        m * v_vec + q * A
    end

    return vcat(x_vec, p_can)
end

function canonical_to_velocity(
        x::AbstractVector, p_can::AbstractVector, p;
        relativistic::Bool = false, t = 0.0
    )
    (; q, m, c, pf) = p
    x_vec = SA[x[1], x[2], x[3]]
    p_vec = SA[p_can[1], p_can[2], p_can[3]]
    A = _eval_field(eltype(x_vec), pf.A, x_vec, t)
    p_kin = p_vec - q * A

    if relativistic
        γ = √(1 + sum(p_kin .^ 2) / (m^2 * c^2))
        return p_kin / (γ * m)
    else
        return p_kin / m
    end
end

function canonical_to_velocity(u::AbstractVector, p; relativistic::Bool = false, t = 0.0)
    return canonical_to_velocity(u[1:3], u[4:6], p; relativistic, t)
end

function canonical_hamiltonian(u::AbstractVector, p, t = 0.0; relativistic::Bool = false)
    (; q, m, c, pf) = p
    x = SA[u[1], u[2], u[3]]
    p_can = SA[u[4], u[5], u[6]]

    A = _eval_field(eltype(x), pf.A, x, t)
    phi = pf.phi isa ZeroField ? zero(eltype(x)) : pf.phi(x, t)
    p_kin = p_can - q * A

    if relativistic
        return √(m^2 * c^4 + c^2 * sum(p_kin .^ 2)) + q * phi
    else
        return sum(p_kin .^ 2) / (2 * m) + q * phi
    end
end
