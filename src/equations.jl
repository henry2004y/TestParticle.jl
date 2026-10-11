# Tracing equations.

get_x(u) = @inbounds SA[u[1], u[2], u[3]]
get_v(u) = @inbounds SA[u[4], u[5], u[6]]

function get_dv(v::AbstractVector{T}, x, p, t) where {T <: AbstractFloat}
    q2m, m, Efunc, Bfunc, Ffunc = p
    E = SVector{3, T}(Efunc(x, t))
    B = SVector{3, T}(Bfunc(x, t))
    F = SVector{3, T}(Ffunc(x, t))

    return T(q2m) * (v × B + E) + F / T(m)
end

function get_dv(v, x, p, t)
    q2m, m, Efunc, Bfunc, Ffunc = p
    E = Efunc(x, t)
    B = Bfunc(x, t)
    F = Ffunc(x, t)

    return q2m * (v × B + E) + F / m
end

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
    q2m, m, Efunc, Bfunc, Ffunc = p
    E = Efunc(x, t)
    B = Bfunc(x, t)
    F = Ffunc(x, t)

    val = q2m * (v × B + E) + F / m
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

"""
    trace_relativistic!(dy, y, p, t)

ODE equations for relativistic charged particle (x, γv) moving in EM field with in-place form.
"""
function trace_relativistic!(dy, y, p, t)
    γv = get_v(y)
    v = get_relativistic_v(γv)
    @inbounds dy[1:3] = v
    @inbounds dy[4:6] = get_dv(v, y, p, t)

    return
end

"""
    trace_gc_exb!(dx, x, p, t)

Equations for tracing the guiding center using the ExB drift and parallel velocity from a reference trajectory.
"""
function trace_gc_exb!(dx, x, p, t)
    # TODO: support external forces
    q2m, _, Efunc, Bfunc, _, sol = p
    xu = sol(t)
    v = get_v(xu)
    E = Efunc(x, t)
    B = Bfunc(x, t)

    Bmag = norm(B)
    b = B / Bmag
    v_par = (v ⋅ b) .* b

    @inbounds dx[1:3] = (E × b) / Bmag + v_par

    return
end

"""
    trace_gc_flr!(dx, x, p, t)

Equations for tracing the guiding center using the ExB drift with FLR corrections and parallel velocity.
"""
function trace_gc_flr!(dx, x, p, t)
    # TODO: support external forces
    q2m, _, Efunc, Bfunc, _, sol = p
    xu = sol(t)
    xp = get_x(xu)
    v = get_v(xu)
    E = Efunc(x, t)
    B = Bfunc(x, t)

    # B at particle position
    Bx = Bfunc(xp, t)
    Bmag_particle = norm(Bx)
    b_particle = Bx / Bmag_particle
    v_par = (v ⋅ b_particle) .* b_particle
    v_perp = v - v_par

    r4 = (v_perp ⋅ v_perp) / (4 * (q2m * Bmag_particle)^2)

    # Helper for FLR term: (E × B) / B²
    EB(x_in) = begin
        E_in = Efunc(x_in, t)
        B_in = Bfunc(x_in, t)
        (E_in × B_in) / (B_in ⋅ B_in)
    end

    # dx = EB(x) + r^2/4 * ∇²(EB) + v_par
    # EB(x) is redundant, use E and B directly
    @inbounds dx[1:3] =
        (E × B) / (B ⋅ B) + r4 * Tensors.laplace.(EB, Tensors.Vec(x...)) + v_par

    return
end

"""
    trace_relativistic(y, p, t) -> SVector{6}

ODE equations for relativistic charged particle (x, γv) moving in static EM field with out-of-place form.
"""
function trace_relativistic(y, p, t)
    γv = get_v(y)
    v = get_relativistic_v(γv)
    dv = get_dv(v, y, p, t)

    return vcat(v, dv)
end

"""
    trace_normalized!(dy, y, p, t)

Normalized ODE equations for charged particle moving in EM field with in-place form.
If the field is in 2D X-Y plane, periodic boundary should be applied for the field in z via
the extrapolation function.
"""
function trace_normalized!(dy, y, p, t)
    v = get_v(y)
    E = get_EField(p)(y, t)
    B = get_BField(p)(y, t)

    @inbounds dy[1:3] = v
    @inbounds dy[4:6] = v × B + E

    return
end

"""
    trace_normalized(y, p, t)

Normalized ODE equations for charged particle moving in EM field with out-of-place form.
"""
function trace_normalized(y, p, t)
    v = get_v(y)
    E = get_EField(p)(y, t)
    B = get_BField(p)(y, t)
    dv = v × B + E

    return vcat(v, dv)
end

"""
    trace_relativistic_normalized!(dy, y, p, t)

Normalized ODE equations for relativistic charged particle (x, γv) moving in EM field with in-place form.
"""
function trace_relativistic_normalized!(dy, y, p, t)
    E = get_EField(p)(y, t)
    B = get_BField(p)(y, t)
    γv = get_v(y)

    v = get_relativistic_v(γv; c = 1)
    @inbounds dy[1:3] = v
    @inbounds dy[4:6] = v × B + E

    return
end

"""
    trace_relativistic_normalized(y, p, t)

Normalized ODE equations for relativistic charged particle (x, γv) moving in EM field with out-of-place form.
"""
function trace_relativistic_normalized(y, p, t)
    E = get_EField(p)(y, t)
    B = get_BField(p)(y, t)
    γv = get_v(y)

    v = get_relativistic_v(γv; c = 1)
    dv = v × B + E

    return vcat(v, dv)
end

"""
    get_B_parameters(x, t, Bfunc)

Evaluate the magnetic field `Bfunc` at `(x, t)` and return the field vector `B`,
its magnitude `Bmag`, the unit vector `b̂ = B / Bmag`, the gradient `∇B`, and the
Jacobian `JB` of the field.
"""
@inline function get_B_parameters(x, t, Bfunc)
    B, JB = _get_B_jacobian(x, t, Bfunc)
    Bmag = norm(B)
    b̂ = B / Bmag
    # Grad-B from Jacobian
    ∇B = JB' * b̂

    return B, Bmag, b̂, ∇B, JB
end

@inline function get_E_parameters(x, t, Efunc)
    E = Efunc(x, t)
    JE = jacobian(Efunc, x, t)

    return E, JE
end

"""
    trace_gc_drifts!(dx, x, p, t)

Equations for tracing the guiding center using analytical drifts, including the grad-B drift, curvature drift, and ExB drift.
Parallel velocity is also added. This expression requires the full particle trajectory `p.sol`.
"""
function trace_gc_drifts!(dx, x, p, t)
    # TODO: support external forces
    q2m, _, Efunc, Bfunc, _, sol = p
    xu = sol(t)
    v = get_v(xu)
    E = Efunc(x, t)

    B, Bmag, b, ∇B, JB = get_B_parameters(x, t, Bfunc)

    v_par_val = v ⋅ b
    v_par = v_par_val .* b
    v_perp = v - v_par
    Ω = q2m * Bmag

    # Curvature vector κ = (b̂ ⋅ ∇) b̂
    # κ = (JB * b̂ - b̂ * (∇B ⋅ b̂)) / Bmag
    κ = (JB * b + b * (-∇B ⋅ b)) / Bmag

    v_E = (E × b) / Bmag
    w = v_perp - v_E

    # w^2*(b×∇|B|)/(2*Ω*B) + v∥^2*(b×κ)/Ω + v_E + v∥
    @inbounds dx[1:3] = (w[1]^2 + w[2]^2 + w[3]^2) * (b × ∇B) / (2 * Ω * Bmag) +
        v_par_val^2 * (b × κ) / Ω +
        v_E + v_par

    return
end


"""
    trace_gc!(dy, y, p, t)

Guiding center equations for nonrelativistic charged particle moving in EM field with in-place form.
Variable `y = (x, y, z, u)`, where `u` is the velocity along the magnetic field at (x,y,z).
"""
function trace_gc!(dy, y, p, t)
    v1, v2, v3, du = get_gc_derivatives(y, p, t)

    @inbounds dy[1] = v1
    @inbounds dy[2] = v2
    @inbounds dy[3] = v3
    @inbounds dy[4] = du

    return
end

"""
    trace_gc(y, p, t)

Guiding center equations for nonrelativistic charged particle moving in EM field with out-of-place form.
"""
function trace_gc(y, p, t)
    v1, v2, v3, du = get_gc_derivatives(y, p, t)
    return SVector{4}(v1, v2, v3, du)
end

"""
    get_gc_velocity(y, p, t)

Get the guiding center velocity.
"""
function get_gc_velocity(y, p, t)
    v1, v2, v3, _ = get_gc_derivatives(y, p, t)
    return SVector{3}(v1, v2, v3)
end

function get_gc_derivatives(y, p, t)
    # TODO: support external forces
    q, q2m, μ, Efunc, Bfunc = p
    X = get_x(y)

    E = Efunc(X, t)
    B, Bmag, b̂, ∇B, JB = get_B_parameters(X, t, Bfunc)

    # ∇ × b̂ = (∇ × B + b̂ × ∇B) / B
    # ∇ × B from JB (Jacobian of B)
    curlB = SVector{3}(JB[3, 2] - JB[2, 3], JB[1, 3] - JB[3, 1], JB[2, 1] - JB[1, 2])
    curlb = (curlB + b̂ × ∇B) / Bmag

    # effective EM fields
    # B* = B + (m/q) u (∇ × b)
    # E* = E - (μ/q) ∇B
    # In CGS: B* = B + (c p_par / q) (∇ × b). In SI, c -> 1.
    Bᵉ = B + (y[4] / q2m) * curlb
    Eᵉ = E - (μ / q) * ∇B

    inv_Bparᵉ = inv(b̂ ⋅ Bᵉ)

    # dx/dt = (p_par/m * B* +  E* × b ) / B*_par
    #       = (u * B* + E* × b) / B*_par
    # In CGS: c/q * q E* x b. In SI, c=1.
    v = (y[4] * Bᵉ + Eᵉ × b̂) * inv_Bparᵉ

    # dp_par/dt = q/B*_par * B* ⋅ E* => du/dt = (q/m)/B*_par * B* ⋅ E*
    du = q2m * inv_Bparᵉ * Bᵉ ⋅ Eᵉ

    return v[1], v[2], v[3], du
end

"""
    trace_fieldline!(dx, x, p, s)

Equation for tracing magnetic field lines with in-place form.
The parameter `p` is the magnetic field function.
Note that the independent variable `s` represents the arc length.
"""
function trace_fieldline!(dx, x, p, s)
    B = p(x, s)
    val = normalize(B)
    @inbounds for i in eachindex(dx, val)
        dx[i] = val[i]
    end
    return
end

"""
    trace_fieldline(x, p, s)

Equation for tracing magnetic field lines with out-of-place form.
"""
function trace_fieldline(x, p, s)
    B = p(x, s)
    return normalize(B)
end

"""
    get_work_rates(xu, p, t)

Calculate the work rates done by the electric field and the betatron acceleration.
Returns a tuple `(P_par, P_fermi, P_grad, P_betatron)`.
"""
@inline function get_work_rates(xu, p, t, magnetic_properties = nothing, E_field = nothing)
    q2m, m, Efunc, Bfunc, _ = p
    r = get_x(xu)
    q = q2m * m

    if E_field === nothing
        E = Efunc(r, t)
    else
        E = E_field
    end

    if magnetic_properties === nothing
        B, ∇B, κ, b̂, Bmag = get_magnetic_properties(r, t, Bfunc)
    else
        B, ∇B, κ, b̂, Bmag = magnetic_properties
    end

    if Bmag == 0
        return SVector{4, eltype(xu)}(0, 0, 0, 0)
    end

    v = get_v(xu)

    # Parallel velocity
    v_par_val = v ⋅ b̂
    v_par = v_par_val .* b̂
    v_perp = v - v_par

    # Magnetic moment
    w_sq = v_perp ⋅ v_perp
    μ = m * w_sq / (2 * Bmag)

    # 1. Parallel Work: q v_par (E ⋅ b)
    P_par = q * v_par_val * (E ⋅ b̂)

    # 2. Fermi Work: m v_par^2 / B (b × κ) ⋅ E
    P_fermi = (m * v_par_val^2 / Bmag) * ((b̂ × κ) ⋅ E)

    # 3. Gradient Drift Work: μ / B (b × ∇B) ⋅ E
    P_grad = (μ / Bmag) * ((b̂ × ∇B) ⋅ E)

    # 4. Betatron Work: μ ∂B/∂t
    dBdt_val = derivative_t(Bfunc, r, t) ⋅ b̂

    P_betatron = μ * dBdt_val

    return SVector{4}(P_par, P_fermi, P_grad, P_betatron)
end

"""
    get_work_rates_gc(xv, p, t)

Calculate the work rates done by the electric field and the betatron acceleration for guiding center.
"""
function get_work_rates_gc(xv, p, t)
    # p = (q, q2m, μ, Efunc, Bfunc)
    q, q2m, μ, Efunc, Bfunc = p
    r = get_x(xv)
    v_par = xv[4]
    E = Efunc(r, t)

    B, ∇B, κ, b̂, Bmag = get_magnetic_properties(r, t, Bfunc)

    if Bmag == 0
        return SVector{4, eltype(xv)}(0, 0, 0, 0)
    end

    m = q / q2m

    # 1. Parallel Work: q v_par (E ⋅ b)
    P_par = q * v_par * (E ⋅ b̂)

    # 2. Fermi Work: m v_par^2 / B (b × κ) ⋅ E
    P_fermi = (m * v_par^2 / Bmag) * ((b̂ × κ) ⋅ E)

    # 3. Gradient Drift Work: μ / B (b × ∇B) ⋅ E
    P_grad = (μ / Bmag) * ((b̂ × ∇B) ⋅ E)

    # 4. Betatron Work: μ ∂B/∂t
    dBdt_val = derivative_t(Bfunc, r, t) ⋅ b̂

    P_betatron = μ * dBdt_val

    return SVector{4}(P_par, P_fermi, P_grad, P_betatron)
end

function trace_canonical!(dy, y, p, t)
    q, m, _, pf = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = pf.A isa ZeroField ? zero(x) : SVector{3}(pf.A(x, t))
    v = (p_can - q * A) / m

    grad_phi = pf.grad_phi isa ZeroField ? zero(x) : SVector{3}(pf.grad_phi(x, t))
    grad_A = pf.grad_A(x, t)

    dp = q * (grad_A * v) - q * grad_phi

    @inbounds dy[1] = v[1]
    @inbounds dy[2] = v[2]
    @inbounds dy[3] = v[3]
    @inbounds dy[4] = dp[1]
    @inbounds dy[5] = dp[2]
    @inbounds dy[6] = dp[3]
    return
end

function trace_canonical(y, p, t)
    q, m, _, pf = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = pf.A isa ZeroField ? zero(x) : SVector{3}(pf.A(x, t))
    v = (p_can - q * A) / m

    grad_phi = pf.grad_phi isa ZeroField ? zero(x) : SVector{3}(pf.grad_phi(x, t))
    grad_A = pf.grad_A(x, t)

    dp = q * (grad_A * v) - q * grad_phi
    return vcat(v, dp)
end

function trace_canonical_relativistic!(dy, y, p, t)
    q, m, c_val, pf = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = pf.A isa ZeroField ? zero(x) : SVector{3}(pf.A(x, t))
    p_kin = p_can - q * A
    γ = √(1 + sum(p_kin .^ 2) / (m^2 * c_val^2))
    v = p_kin / (γ * m)

    grad_phi = pf.grad_phi isa ZeroField ? zero(x) : SVector{3}(pf.grad_phi(x, t))
    grad_A = pf.grad_A(x, t)

    dp = q * (grad_A * v) - q * grad_phi

    @inbounds dy[1] = v[1]
    @inbounds dy[2] = v[2]
    @inbounds dy[3] = v[3]
    @inbounds dy[4] = dp[1]
    @inbounds dy[5] = dp[2]
    @inbounds dy[6] = dp[3]
    return
end

function trace_canonical_relativistic(y, p, t)
    q, m, c_val, pf = p
    x = @inbounds SA[y[1], y[2], y[3]]
    p_can = @inbounds SA[y[4], y[5], y[6]]

    A = pf.A isa ZeroField ? zero(x) : SVector{3}(pf.A(x, t))
    p_kin = p_can - q * A
    γ = √(1 + sum(p_kin .^ 2) / (m^2 * c_val^2))
    v = p_kin / (γ * m)

    grad_phi = pf.grad_phi isa ZeroField ? zero(x) : SVector{3}(pf.grad_phi(x, t))
    grad_A = pf.grad_A(x, t)

    dp = q * (grad_A * v) - q * grad_phi
    return vcat(v, dp)
end

function velocity_to_canonical(
        x0::AbstractVector, v0::AbstractVector, p;
        relativistic::Bool = false, t = 0.0
    )
    q, m, c_val, pf = p
    x_vec = SA[x0[1], x0[2], x0[3]]
    v_vec = SA[v0[1], v0[2], v0[3]]
    A = pf.A isa ZeroField ? zero(x_vec) : SVector{3}(pf.A(x_vec, t))

    p_can = if relativistic
        v2 = sum(v_vec .^ 2)
        γ = 1 / √(1 - v2 / c_val^2)
        γ * m * v_vec + q * A
    else
        m * v_vec + q * A
    end

    if x0 isa StaticVector && v0 isa StaticVector
        return vcat(x_vec, p_can)
    else
        return [x_vec..., p_can...]
    end
end

function canonical_to_velocity(
        x::AbstractVector, p_can::AbstractVector, p;
        relativistic::Bool = false, t = 0.0
    )
    q, m, c_val, pf = p
    x_vec = SA[x[1], x[2], x[3]]
    p_vec = SA[p_can[1], p_can[2], p_can[3]]
    A = pf.A isa ZeroField ? zero(x_vec) : SVector{3}(pf.A(x_vec, t))
    p_kin = p_vec - q * A

    if relativistic
        γ = √(1 + sum(p_kin .^ 2) / (m^2 * c_val^2))
        return p_kin / (γ * m)
    else
        return p_kin / m
    end
end

function canonical_to_velocity(u::AbstractVector, p; relativistic::Bool = false, t = 0.0)
    return canonical_to_velocity(u[1:3], u[4:6], p; relativistic, t)
end

function canonical_hamiltonian(u::AbstractVector, p, t = 0.0; relativistic::Bool = false)
    q, m, c_val, pf = p
    x = SA[u[1], u[2], u[3]]
    p_can = SA[u[4], u[5], u[6]]

    A = pf.A isa ZeroField ? zero(x) : SVector{3}(pf.A(x, t))
    phi = pf.phi isa ZeroField ? zero(eltype(x)) : pf.phi(x, t)
    p_kin = p_can - q * A

    if relativistic
        return √(m^2 * c_val^4 + c_val^2 * sum(p_kin .^ 2)) + q * phi
    else
        return sum(p_kin .^ 2) / (2 * m) + q * phi
    end
end

