# Guiding-center equations.

"""
    trace_gc_exb!(dx, x, p, t)

Equations for tracing the guiding center using the ExB drift and parallel velocity from a reference trajectory.
"""
function trace_gc_exb!(dx, x, p, t)
    # TODO: support external forces
    v = get_v(p.sol(t))
    E = p.E(x, t)
    B = p.B(x, t)

    Bmag = norm(B)
    b = B / Bmag

    @inbounds dx[1:3] = (E × b) / Bmag + (v ⋅ b) * b

    return
end

"""
    trace_gc_flr!(dx, x, p, t)

Equations for tracing the guiding center using the ExB drift with FLR corrections and parallel velocity.
"""
function trace_gc_flr!(dx, x, p, t)
    # TODO: support external forces
    (; q2m, sol) = p
    Efunc, Bfunc = p.E, p.B
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
        (E × B) / (B ⋅ B) + r4 * vector_laplacian(EB, SVector{3}(x)) + v_par

    return
end

# Row 3(j - 1) + i, column k of the nested Jacobian is ∂²fᵢ/∂xⱼ∂xₖ.
function vector_laplacian(f, x::SVector{3})
    H = ForwardDiff.jacobian(y -> vec(ForwardDiff.jacobian(f, y)), x)
    return SVector(ntuple(i -> H[i, 1] + H[3 + i, 2] + H[6 + i, 3], 3))
end

"""
    trace_gc_drifts!(dx, x, p, t)

Equations for tracing the guiding center using analytical drifts, including the grad-B drift, curvature drift, and ExB drift.
Parallel velocity is also added. This expression requires the full particle trajectory `p.sol`.
"""
function trace_gc_drifts!(dx, x, p, t)
    # TODO: support external forces
    (; q2m, sol) = p
    v = get_v(sol(t))
    E = p.E(x, t)
    (; ∇B, κ, b̂, Bmag) = get_magnetic_properties(x, t, p.B)

    v_par_val = v ⋅ b̂
    v_par = v_par_val * b̂
    Ω = q2m * Bmag
    v_E = (E × b̂) / Bmag
    w = v - v_par - v_E

    # w^2*(b×∇|B|)/(2*Ω*B) + v∥^2*(b×κ)/Ω + v_E + v∥
    @inbounds dx[1:3] = (w ⋅ w) * (b̂ × ∇B) / (2 * Ω * Bmag) +
        v_par_val^2 * (b̂ × κ) / Ω +
        v_E + v_par

    return
end

@inline _gc_gamma(::GuidingCenter, u, Bmag, p) = one(u)

@inline function _gc_gamma(::RelativisticGuidingCenter, u, Bmag, p)
    u² = u^2 + 2 * p.μ * Bmag * (p.q2m / p.q)
    c_val = oftype(float(u²), _speed_of_light(p))
    return sqrt(1 + u² / c_val^2)
end

@muladd @inline function _trace_gc(y, model, p, t)
    # TODO: support external forces
    (; q, q2m, μ) = p
    X = get_x(y)

    E = p.E(X, t)
    B, JB = _get_B_jacobian(X, t, p.B)
    Bmag = norm(B)
    if iszero(Bmag)
        return zero(SVector{4, eltype(B)})
    end
    inv_Bmag = inv(Bmag)
    b̂ = B * inv_Bmag
    ∇B = JB' * b̂

    # ∇ × b̂ = (∇ × B + b̂ × ∇B) / B
    # ∇ × B from JB (Jacobian of B)
    curlB = SVector{3}(JB[3, 2] - JB[2, 3], JB[1, 3] - JB[3, 1], JB[2, 1] - JB[1, 2])
    curlb = (curlB + b̂ × ∇B) * inv_Bmag

    # Effective EM fields:
    # B* = B + (m/q) u (∇ × b̂)
    # E* = E - (μ / (q γ)) ∇B
    # where u = v∥ (non-relativistic, γ = 1) or u = γv∥ (relativistic).
    u = y[4]
    inv_γ = inv(_gc_gamma(model, u, Bmag, p))
    Bᵉ = B + (u / q2m) * curlb
    Eᵉ = E - (μ * inv_γ / q) * ∇B

    inv_Bparᵉ = inv(b̂ ⋅ Bᵉ)

    # dX/dt = ((u / γ) B* + E* × b̂) / B*_par
    v = ((u * inv_γ) * Bᵉ + Eᵉ × b̂) * inv_Bparᵉ

    # du/dt = (q/m) (B* ⋅ E*) / B*_par
    du = q2m * inv_Bparᵉ * (Bᵉ ⋅ Eᵉ)

    return SVector{4}(v[1], v[2], v[3], du)
end

"""
    trace_gc!(dy, y, p, t)

Guiding center equations for nonrelativistic charged particle moving in EM field with in-place form.
Variable `y = (x, y, z, u)`, where `u` is the velocity along the magnetic field at (x,y,z).
"""
trace_gc!(dy, y, p, t) = _store!(dy, @inline trace_gc(y, p, t))

"""
    trace_gc(y, p, t)

Guiding center equations for nonrelativistic charged particle moving in EM field with out-of-place form.
"""
trace_gc(y, p, t) = _trace_gc(y, GuidingCenter(), p, t)

"""
    trace_gc_relativistic!(dy, y, p, t)

Relativistic guiding-center equations in in-place form for state `y = (X_x, X_y, X_z, u_∥)`,
where `u_∥ = γ v_∥` is the parallel proper velocity and `p.μ = p_⊥² / (2 m B)` is the
relativistic magnetic moment. The speed of light is `p.c` when present, `c` in SI otherwise.
"""
trace_gc_relativistic!(dy, y, p, t) =
    _store!(dy, @inline trace_gc_relativistic(y, p, t))

"""
    trace_gc_relativistic(y, p, t)

Relativistic guiding-center equations in out-of-place form for state `y = (X_x, X_y, X_z, u_∥)`.
"""
trace_gc_relativistic(y, p, t) = _trace_gc(y, RelativisticGuidingCenter(), p, t)

"""
    get_gc_velocity(y, p, t; model = GuidingCenter())

Get the guiding center velocity for `model` (`GuidingCenter()` or
`RelativisticGuidingCenter()`).
"""
get_gc_velocity(y, p, t; model = GuidingCenter()) = get_x(_trace_gc(y, model, p, t))
