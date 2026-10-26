# Guiding-center equations.

"""
    trace_gc_exb!(dx, x, p, t)

Equations for tracing the guiding center using the ExB drift and parallel velocity from a reference trajectory.
"""
function trace_gc_exb!(dx, x, p, t)
    # TODO: support external forces
    (; q2m, sol) = p
    Efunc, Bfunc = p.E, p.B
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
    Efunc, Bfunc = p.E, p.B
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
    (; q, q2m, μ) = p
    X = get_x(y)

    E = p.E(X, t)
    B, Bmag, b̂, ∇B, JB = get_B_parameters(X, t, p.B)

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
