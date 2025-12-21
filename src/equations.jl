# Tracing equations.

get_x(u) = @inbounds SA[u[1], u[2], u[3]]
get_v(u) = @inbounds SA[u[4], u[5], u[6]]

function get_fields(x, t, Efunc, Bfunc)
   E = Efunc(x, t)
   B = Bfunc(x, t)
   return E, B
end

function get_B_properties(x, t, Bfunc)
   B = Bfunc(x, t)
   Bmag = norm(B)
   b̂ = B / Bmag
   # Gradient of B magnitude
   ∇B = ForwardDiff.gradient(x -> norm(Bfunc(x, t)), x)
   return B, Bmag, b̂, ∇B
end

function get_b_jacobian(x, t, Bfunc)
   bfunc(x) = normalize(Bfunc(x, t))
   ∇b̂ = ForwardDiff.jacobian(bfunc, x)
   return ∇b̂
end

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

   return nothing
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
   E, B = get_fields(x, t, Efunc, Bfunc)

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
   E, B = get_fields(x, t, Efunc, Bfunc)

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
   dx[1:3] = EB(x) + r4 * Tensors.laplace.(EB, Tensors.Vec(x...)) + v_par

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

"""
     trace_gc_drifts!(dx, x, p, t)

Equations for tracing the guiding center using analytical drifts, including the grad-B drift, curvature drift, and ExB drift.
Parallel velocity is also added. This expression requires the full particle trajectory `p.sol`.
"""
function trace_gc_drifts!(dx, x, p, t)
   q2m, _, Efunc, Bfunc, _, sol = p
   xu = sol(t)
   v = get_v(xu)
   E, B = get_fields(x, t, Efunc, Bfunc)

   # Gradient of B magnitude and other properties
   _, Bmag, b, ∇B = get_B_properties(x, t, Bfunc)

   v_par = (v ⋅ b) .* b
   v_perp = v - v_par
   Ω = q2m * Bmag

   ∇b̂ = get_b_jacobian(x, t, Bfunc)
   κ = ∇b̂ * b  # curvature vector (b.∇)b

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
   q, m, μ, Efunc, Bfunc = p
   q2m = q / m
   X = get_x(y)

   E, B = get_fields(X, t, Efunc, Bfunc)
   _, _, b̂, ∇B = get_B_properties(X, t, Bfunc)
   ∇b̂ = get_b_jacobian(X, t, Bfunc)

   function E2B(x)
      E² = Efunc(x, t) ⋅ Efunc(x, t)
      B² = Bfunc(x, t) ⋅ Bfunc(x, t)
      vE² = E² / B²
   end

   ∇vE² = SVector{3}(ForwardDiff.gradient(E2B, X))

   # ∇ × b̂
   curlb = SVector{3}(∇b̂[3, 2] - ∇b̂[2, 3], ∇b̂[1, 3] - ∇b̂[3, 1], ∇b̂[2, 1] - ∇b̂[1, 2])
   # effective EM fields
   Eᵉ = @. E - (μ * ∇B - 0.5 * m * ∇vE²) / q
   Bᵉ = @. B + y[4] / q2m * curlb
   Bparᵉ = b̂ ⋅ Bᵉ # effective B field parallel to B

   dy[1] = (y[4] * Bᵉ[1] + Eᵉ[2] * b̂[3] - Eᵉ[3] * b̂[2]) / Bparᵉ
   dy[2] = (y[4] * Bᵉ[2] + Eᵉ[3] * b̂[1] - Eᵉ[1] * b̂[3]) / Bparᵉ
   dy[3] = (y[4] * Bᵉ[3] + Eᵉ[1] * b̂[2] - Eᵉ[2] * b̂[1]) / Bparᵉ
   dy[4] = q2m / Bparᵉ * Bᵉ ⋅ Eᵉ

   return
end

"""
1st order approximation of guiding center equations.
"""
function trace_gc_1st!(dy, y, p::GCTuple, t)
   q, m, μ, Efunc, Bfunc = p
   q2m = q / m
   X = get_x(y)

   E, B = get_fields(X, t, Efunc, Bfunc)
   _, Bmag, b̂, ∇B = get_B_properties(X, t, Bfunc)
   ∇b̂ = get_b_jacobian(X, t, Bfunc)

   Ω = q * Bmag / m
   u = y[4]

   # effective EM fields
   Eᵉ = @. E - (μ * ∇B) / q
   κ = ∇b̂ * b̂  # curvature
   vX = u * b̂ + b̂ × (κ * u^2 - q2m * Eᵉ) / Ω

   dy[1] = vX[1]
   dy[2] = vX[2]
   dy[3] = vX[3]
   dy[4] = q2m * (b̂ + u * (b̂ × κ) / Ω) ⋅ Eᵉ

   return
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
