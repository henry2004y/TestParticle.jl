# Tracing equations.

get_x(u) = @inbounds SA[u[1], u[2], u[3]]
get_v(u) = @inbounds SA[u[4], u[5], u[6]]

function get_dv(x, v, p, t)
   q2m, m, Efunc, Bfunc, Ffunc = p
   E = Efunc(x, t)
   B = Bfunc(x, t)
   F = Ffunc(x, t)

   SVector{3}(q2m * (v × B + E) + F / m)
end

function get_relativistic_v(γv; c = c)
   γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
   if γ²v² > eps(eltype(γv))
      v̂ = normalize(γv)
   else # no velocity
      v̂ = SVector{3, eltype(γv)}(0, 0, 0)
   end
   √(γ²v² / (1 + γ²v² / c^2)) * v̂
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
   vcat(v, dv)
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

   vcat(v, dv)
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

   vcat(v, dv)
end

@inline function get_B_parameters(x, t, Bfunc)
   # Compute B and its Jacobian in a single pass using ForwardDiff
   result = DiffResults.JacobianResult(x)
   result = ForwardDiff.jacobian!(result, x -> Bfunc(x, t), x)

   B = SVector{3}(DiffResults.value(result))
   JB = SMatrix{3, 3}(DiffResults.jacobian(result))

   Bmag = norm(B)
   b̂ = B / Bmag

   # ∇|B| = (J_B' * b̂)
   ∇B = JB' * b̂

   return B, Bmag, b̂, ∇B, JB
end

@inline function get_E_parameters(x, t, Efunc)
   result = DiffResults.JacobianResult(x)
   result = ForwardDiff.jacobian!(result, x -> Efunc(x, t), x)

   E = SVector{3}(DiffResults.value(result))
   JE = SMatrix{3, 3}(DiffResults.jacobian(result))

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

   # v⟂^2*(b×∇|B|)/(2*Ω*B) + v∥^2*(b×κ)/Ω + (E×b)/B + v∥
   @inbounds dx[1:3] = norm(v_perp)^2 * (b × ∇B) / (2 * Ω * Bmag) +
                       norm(v_par)^2 * (b × κ) / Ω +
                       (E × b) / Bmag + v_par

   return
end

"""
    trace_gc!(dy, y, p, t)

Guiding center equations for nonrelativistic charged particle moving in EM field with in-place form.
Variable `y = (x, y, z, u)`, where `u` is the velocity along the magnetic field at (x,y,z).
"""
function trace_gc!(dy, y, p::GCTuple, t)
   v, du = get_gc_derivatives(y, p, t)

   @inbounds dy[1:3] = v
   @inbounds dy[4] = du

   return
end

"""
    get_gc_velocity(y, p, t)

Get the guiding center velocity vector for state `y`, parameters `p` at time `t`.
"""
function get_gc_velocity(y, p::GCTuple, t)
   v, _ = get_gc_derivatives(y, p, t)
   return v
end

function get_gc_derivatives(y, p::GCTuple, t)
   q, m, μ, Efunc, Bfunc = p
   q2m = q / m
   X = get_x(y)

   E, JE = get_E_parameters(X, t, Efunc)
   B, Bmag, b̂, ∇B, JB = get_B_parameters(X, t, Bfunc)

   # ∇ × b̂ = (∇ × B + b̂ × ∇B) / B
   # ∇ × B from JB
   curlB = SVector{3}(JB[3, 2] - JB[2, 3], JB[1, 3] - JB[3, 1], JB[2, 1] - JB[1, 2])
   curlb = (curlB + b̂ × ∇B) / Bmag

   # ∇(E²/B²) = 2/B² * (JE'*E - (E²/B)*∇B)
   # Uses JE (Jacobian of E) and JB (via ∇B)
   E² = E ⋅ E
   ∇vE² = (2 / Bmag^2) * (JE' * E - (E² / Bmag) * ∇B)

   # effective EM fields
   Eᵉ = @. E - (μ * ∇B - 0.5 * m * ∇vE²) / q
   Bᵉ = @. B + y[4] / q2m * curlb
   Bparᵉ = b̂ ⋅ Bᵉ # effective B field parallel to B

   v1 = (y[4] * Bᵉ[1] + Eᵉ[2] * b̂[3] - Eᵉ[3] * b̂[2]) / Bparᵉ
   v2 = (y[4] * Bᵉ[2] + Eᵉ[3] * b̂[1] - Eᵉ[1] * b̂[3]) / Bparᵉ
   v3 = (y[4] * Bᵉ[3] + Eᵉ[1] * b̂[2] - Eᵉ[2] * b̂[1]) / Bparᵉ
   du = q2m / Bparᵉ * Bᵉ ⋅ Eᵉ

   SVector{3}(v1, v2, v3), du
end

"""
1st order approximation of guiding center equations.
"""
function trace_gc_1st!(dy, y, p::GCTuple, t)
   v, du = get_gc_1st_derivatives(y, p, t)

   @inbounds dy[1:3] = v
   @inbounds dy[4] = du

   return
end

"""
    get_gc_1st_velocity(y, p, t)

Get the 1st order guiding center velocity for state `y`, parameters `p` at time `t`.
"""
function get_gc_1st_velocity(y, p::GCTuple, t)
   v, _ = get_gc_1st_derivatives(y, p, t)
   return v
end

function get_gc_1st_derivatives(y, p::GCTuple, t)
   q, m, μ, Efunc, Bfunc = p
   q2m = q / m
   X = get_x(y)

   E = Efunc(X, t)
   B, Bmag, b̂, ∇B, JB = get_B_parameters(X, t, Bfunc)

   Ω = q * Bmag / m
   u = y[4]

   # effective EM fields
   Eᵉ = @. E - (μ * ∇B) / q

   # Curvature vector κ = (b̂ ⋅ ∇) b̂
   κ = (JB * b̂ + b̂ * (-∇B ⋅ b̂)) / Bmag

   vX = u * b̂ + b̂ × (κ * u^2 - q2m * Eᵉ) / Ω

   du = q2m * (b̂ + u * (b̂ × κ) / Ω) ⋅ Eᵉ

   vX, du
end

"""
    trace_fieldline!(dx, x, p, s)

Equation for tracing magnetic field lines with in-place form.
The parameter `p` is the magnetic field function.
Note that the independent variable `s` represents the arc length.
"""
function trace_fieldline!(dx, x, p, s)
   B = p(x, s)
   dx .= normalize(B)
end

"""
    trace_fieldline(x, p, s)

Equation for tracing magnetic field lines with out-of-place form.
"""
function trace_fieldline(x, p, s)
   B = p(x, s)
   normalize(B)
end
