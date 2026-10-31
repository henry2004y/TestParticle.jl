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

"""
    trace_gc!(dy, y, p, t)

Guiding center equations for nonrelativistic charged particle moving in EM field with in-place form.
Variable `y = (x, y, z, u)`, where `u` is the velocity along the magnetic field at (x,y,z).
"""
trace_gc!(dy, y, p, t) = _store!(dy, @inline trace_gc(y, p, t))

"""
    get_gc_velocity(y, p, t)

Get the guiding center velocity.
"""
get_gc_velocity(y, p, t) = get_x(trace_gc(y, p, t))

"""
    trace_gc(y, p, t)

Guiding center equations for nonrelativistic charged particle moving in EM field with out-of-place form.
"""
function trace_gc(y, p, t)
    return _trace_gc(y, p, t, get_magnetic_properties(get_x(y), t, p.B))
end

"""
    _trace_gc(y, p, t, mp)

Same as [`trace_gc`](@ref) with the [`get_magnetic_properties`](@ref) output
`mp` at `get_x(y)` supplied by the caller, so that a step and the adiabaticity
check can share one evaluation of the field derivatives.
"""
@inline function _trace_gc(y, p, t, mp)
    # TODO: support external forces
    (; q, q2m, μ) = p
    X = get_x(y)

    E = p.E(X, t)
    (; B, Bmag, b̂, ∇B, JB) = mp

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

    return SVector{4}(v[1], v[2], v[3], du)
end
