# Tracing equations.

get_x(u) = @inbounds SA[u[1], u[2], u[3]]
get_v(u) = @inbounds SA[u[4], u[5], u[6]]

function get_dv(x, v, p, t)
    q2m, m, Efunc, Bfunc, Ffunc = p
    E = Efunc(x, t)
    B = Bfunc(x, t)
    F = Ffunc(x, t)

    return SVector{3}(q2m * (v × B + E) + F / m)
end

function get_relativistic_v(γv; c = c)
    γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
    if γ²v² > eps(eltype(γv))
        v̂ = normalize(γv)
    else # no velocity
        v̂ = SVector{3, eltype(γv)}(0, 0, 0)
    end
    return √(γ²v² / (1 + γ²v² / c^2)) * v̂
end

"""
    trace!(dy, y, p, t)

ODE equations for charged particle moving in EM field and external force field with in-place form.
"""
function trace!(dy, y, p, t)
    v = get_v(y)
    @inbounds dy[1:3] = v
    @inbounds dy[4:6] = get_dv(y, v, p, t)

    return
end

"""
    trace(y, p, t)::SVector{6}

ODE equations for charged particle moving in EM field and external force field with out-of-place form.
"""
function trace(y, p, t)
    v = y[SA[4:6...]]
    dv = get_dv(y, v, p, t)
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
    @inbounds dy[4:6] = get_dv(y, v, p, t)

    return
end

"""
    trace_gc_exb!(dx, x, p, t)

Equations for tracing the guiding center using the ExB drift and parallel velocity from a reference trajectory.
"""
function trace_gc_exb!(dx, x, p, t)
    q2m, _, Efunc, Bfunc, _, sol = p
    xu = sol(t)
    v = get_v(xu)
    E = Efunc(x, t)
    B = Bfunc(x, t)

    Bmag = norm(B)
    b = B / Bmag
    v_par = (v ⋅ b) .* b

    dx[1:3] = (E × b) / Bmag + v_par

    return
end

"""
    trace_gc_flr!(dx, x, p, t)

Equations for tracing the guiding center using the ExB drift with FLR corrections and parallel velocity.
"""
function trace_gc_flr!(dx, x, p, t)
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

    r4 = (norm(v_perp) / q2m / Bmag_particle)^2 / 4

    # Helper for FLR term: (E × B) / B²
    EB(x_in) = begin
        E_in = Efunc(x_in, t)
        B_in = Bfunc(x_in, t)
        (E_in × B_in) / (B_in ⋅ B_in)
    end

    # dx = EB(x) + r^2/4 * ∇²(EB) + v_par
    # EB(x) is redundant, use E and B directly
    dx[1:3] = (E × B) / (B ⋅ B) + r4 * Tensors.laplace.(EB, Tensors.Vec(x...)) + v_par

    return
end

"""
    trace_relativistic(y, p, t) -> SVector{6}

ODE equations for relativistic charged particle (x, γv) moving in static EM field with out-of-place form.
"""
function trace_relativistic(y, p, t)
    γv = get_v(y)
    v = get_relativistic_v(γv)
    dv = get_dv(y, v, p, t)

    return vcat(v, dv)
end

"""
    trace_normalized!(dy, y, p, t)

Normalized ODE equations for charged particle moving in EM field with in-place form.
If the field is in 2D X-Y plane, periodic boundary should be applied for the field in z via
the extrapolation function provided by Interpolations.jl.
"""
function trace_normalized!(dy, y, p, t)
    v = get_v(y)
    E = get_EField(p)(y, t)
    B = get_BField(p)(y, t)

    @inbounds dy[1:3] = v
    @inbounds dy[4:6] = SVector{3}(v × B + E)

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
    dv = SVector{3}(v × B + E)

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
    @inbounds dy[4:6] = SVector{3}(v × B + E)

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
    dv = SVector{3}(v × B + E)

    return vcat(v, dv)
end

@inline function get_B_parameters(x, t, Bfunc)
    # Compute B and its Jacobian in a single pass using ForwardDiff
    JB = ForwardDiff.jacobian(r -> Bfunc(r, t), x)
    B = Bfunc(x, t)

    Bmag = norm(B)
    b̂ = B / Bmag

    # ∇|B| = (J_B' * b̂)
    ∇B = JB' * b̂

    return B, Bmag, b̂, ∇B, JB
end

@inline function get_E_parameters(x, t, Efunc)
    JE = ForwardDiff.jacobian(r -> Efunc(r, t), x)
    E = Efunc(x, t)

    return E, JE
end

"""
    trace_gc_drifts!(dx, x, p, t)

Equations for tracing the guiding center using analytical drifts, including the grad-B drift, curvature drift, and ExB drift.
Parallel velocity is also added. This expression requires the full particle trajectory `p.sol`.
"""
function trace_gc_drifts!(dx, x, p, t)
    q2m, _, Efunc, Bfunc, _, sol = p
    xu = sol(t)
    v = get_v(xu)
    E = Efunc(x, t)

    B, Bmag, b, ∇B, JB = get_B_parameters(x, t, Bfunc)

    v_par = (v ⋅ b) .* b
    v_perp = v - v_par
    Ω = q2m * Bmag

    # Curvature vector κ = (b̂ ⋅ ∇) b̂
    # κ = (JB * b̂ - b̂ * (∇B ⋅ b̂)) / Bmag
    κ = (JB * b + b * (-∇B ⋅ b)) / Bmag

    v_E = (E × b) / Bmag
    w = v_perp - v_E

    # w^2*(b×∇|B|)/(2*Ω*B) + v∥^2*(b×κ)/Ω + v_E + v∥
    @inbounds dx[1:3] = norm(w)^2 * (b × ∇B) / (2 * Ω * Bmag) +
        norm(v_par)^2 * (b × κ) / Ω +
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
    return dx .= normalize(B)
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
function get_work_rates(xu, p, t)
    q2m, m, Efunc, Bfunc, _ = p
    r = get_x(xu)
    v = get_v(xu)
    q = q2m * m
    E = Efunc(r, t)

    B, ∇B, κ, b̂, Bmag = get_magnetic_properties(r, t, Bfunc)

    if Bmag == 0
        return SVector{4, eltype(xu)}(0, 0, 0, 0)
    end

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
    dBdt_val = ForwardDiff.derivative(t -> norm(Bfunc(r, t)), t)

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
    dBdt_val = ForwardDiff.derivative(t -> norm(Bfunc(r, t)), t)

    P_betatron = μ * dBdt_val

    return SVector{4}(P_par, P_fermi, P_grad, P_betatron)
end
