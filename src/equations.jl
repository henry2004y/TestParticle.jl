# Tracing equations.

"""
    trace!(dy, y, p::TPTuple, t)
    trace!(dy, y, p::FullTPTuple, t)

ODE equations for charged particle moving in static EM field with in-place form.

ODE equations for charged particle moving in static EM field and external force field with
in-place form.
"""
function trace!(dy, y, p::TPTuple, t)
   q2m, E, B = p

   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   dy[1], dy[2], dy[3] = vx, vy, vz
   # q/m*(E + v × B)
   dy[4] = q2m*(vy*Bz - vz*By + Ex)
   dy[5] = q2m*(vz*Bx - vx*Bz + Ey)
   dy[6] = q2m*(vx*By - vy*Bx + Ez)

   return
end

function trace!(dy, y, p::FullTPTuple, t)
   q, m, E, B, F = p

   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)
   Fx, Fy, Fz = F(y, t)

   dy[1], dy[2], dy[3] = vx, vy, vz
   dy[4] = (q*(vy*Bz - vz*By + Ex) + Fx) / m
   dy[5] = (q*(vz*Bx - vx*Bz + Ey) + Fy) / m
   dy[6] = (q*(vx*By - vy*Bx + Ez) + Fz) / m

   return
end

"""
    trace(y, p::TPTuple, t) -> SVector{6, Float64}
    trace(y, p::FullTPTuple, t) -> SVector{6, Float64}

ODE equations for charged particle moving in static EM field with out-of-place form.

ODE equations for charged particle moving in static EM field and external force field with
out-of-place form.
"""
function trace(y, p::TPTuple, t)
   q2m, E, B = p
   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   dx, dy, dz = vx, vy, vz
   # q/m*(E + v × B)
   dux = q2m*(vy*Bz - vz*By + Ex)
   duy = q2m*(vz*Bx - vx*Bz + Ey)
   duz = q2m*(vx*By - vy*Bx + Ez)
   SVector{6}(dx, dy, dz, dux, duy, duz)
end

function trace(y, p::FullTPTuple, t)
   q, m, E, B, F = p

   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)
   Fx, Fy, Fz = F(y, t)

   dx, dy, dz = vx, vy, vz
   dux = (q*(vy*Bz - vz*By + Ex) + Fx) / m
   duy = (q*(vz*Bx - vx*Bz + Ey) + Fy) / m
   duz = (q*(vx*By - vy*Bx + Ez) + Fz) / m

   SVector{6}(dx, dy, dz, dux, duy, duz)
end

"""
    trace_relativistic!(dy, y, p::TPTuple, t)

ODE equations for relativistic charged particle (x, γv) moving in static EM field with in-place form.
"""
function trace_relativistic!(dy, y, p::TPTuple, t)
   q2m, E, B = p
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   γv = @views SVector{3, eltype(dy)}(y[4:6])
   γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
   if γ²v² > eps(eltype(dy))
      v̂ = normalize(γv)
   else # no velocity
      v̂ = SVector{3, eltype(dy)}(0, 0, 0)
   end
   vmag = √(γ²v² / (1 + γ²v²/c^2))
   vx, vy, vz = vmag * v̂

   dy[1], dy[2], dy[3] = vx, vy, vz
   dy[4] = q2m * (vy*Bz - vz*By + Ex)
   dy[5] = q2m * (vz*Bx - vx*Bz + Ey)
   dy[6] = q2m * (vx*By - vy*Bx + Ez) 

   return
end

"""
    trace_relativistic(y, p::TPTuple, t) -> SVector{6}

ODE equations for relativistic charged particle (x, γv) moving in static EM field with out-of-place form.
"""
function trace_relativistic(y, p::TPTuple, t)
   q2m, E, B = p
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   γv = @views SVector{3, eltype(y)}(y[4:6])
   γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
   if γ²v² > eps(eltype(y))
      v̂ = normalize(γv)
   else # no velocity
      v̂ = SVector{3, eltype(y)}(0, 0, 0)
   end
   vmag = √(γ²v² / (1 + γ²v²/c^2))
   vx, vy, vz = vmag * v̂

   dx, dy, dz = vx, vy, vz
   dux = q2m * (vy*Bz - vz*By + Ex)
   duy = q2m * (vz*Bx - vx*Bz + Ey)
   duz = q2m * (vx*By - vy*Bx + Ez)

   SVector{6}(dx, dy, dz, dux, duy, duz)
end

"""
    trace_normalized!(dy, y, p::TPNormalizedTuple, t)

Normalized ODE equations for charged particle moving in static EM field with in-place form.
If the field is in 2D X-Y plane, periodic boundary should be applied for the field in z via
the extrapolation function provided by Interpolations.jl.
"""
function trace_normalized!(dy, y, p::TPNormalizedTuple, t)
   _, E, B = p

   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   dy[1], dy[2], dy[3] = vx, vy, vz
   # E + v × B
   dy[4] = vy*Bz - vz*By + Ex
   dy[5] = vz*Bx - vx*Bz + Ey
   dy[6] = vx*By - vy*Bx + Ez

   return
end

"""
    trace_relativistic_normalized!(dy, y, p::TPNormalizedTuple, t)

Normalized ODE equations for relativistic charged particle (x, γv) moving in static EM field with in-place form.
"""
function trace_relativistic_normalized!(dy, y, p::TPNormalizedTuple, t)
   _, E, B = p
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   γv = @views SVector{3, eltype(dy)}(y[4:6])
   γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
   if γ²v² > eps(eltype(dy))
      v̂ = normalize(γv)
   else # no velocity
      v̂ = SVector{3, eltype(dy)}(0, 0, 0)
   end
   vmag = √(γ²v² / (1 + γ²v²))
   vx, vy, vz = vmag * v̂

   dy[1], dy[2], dy[3] = vx, vy, vz
   dy[4] = vy*Bz - vz*By + Ex
   dy[5] = vz*Bx - vx*Bz + Ey
   dy[6] = vx*By - vy*Bx + Ez

   return
end

"""
    trace_gc_drifts!(dx, x, p, t)

Equations for tracing the guiding center using analytical drifts, including the grad-B
drift, the curvature drift, the ExB drift. Parallel velocity is also added. This expression
requires the full particle trajectory `p.sol`.
"""
function trace_gc_drifts!(dx, x, p, t)
   q2m, E, B, sol = p
   xu = sol(t)
   v = @views SVector{3, eltype(dx)}(xu[4:6])
   abs_B(x) = norm(B(x))
   gradient_B = ForwardDiff.gradient(abs_B, x)
   Bv = B(x)
   b = normalize(Bv)
   v_par = (v ⋅ b) .* b
   v_perp = v - v_par
   Ω = q2m*norm(Bv)
   κ = ForwardDiff.jacobian(B, x)*Bv  # B⋅∇B
   ## v⟂^2*(B×∇|B|)/(2*Ω*B^2) + v∥^2*(B×(B⋅∇B))/(Ω*B^3) + (E×B)/B^2 + v∥
   dx[1:3] = norm(v_perp)^2*(Bv × gradient_B)/(2*Ω*norm(Bv)^2) +
      norm(v_par)^2*(Bv × κ)/Ω/norm(Bv)^3 + (E(x) × Bv)/norm(Bv)^2 + v_par
end

"""
    trace_gc!(dy, y, p::TPTuple, t)

Guiding center equations for nonrelativistic charged particle moving in static EM field with
in-place form. Variable `y = (x, y, z, u)`, where `u` is the velocity along the magnetic
field at (x,y,z).
"""
function trace_gc!(dy, y, p::GCTuple, t)
   q, m, μ, Efunc, Bfunc = p
   q2m = q / m
   X = @view y[1:3]
   E = Efunc(X, t)
   B = Bfunc(X, t)
   b̂ = normalize(B) # unit B field at X

   Bmag(x) = √(Bfunc(x) ⋅ Bfunc(x))
   ∇B = SVector{3}(ForwardDiff.gradient(Bmag, X))

   function E2B(x)
      E² = Efunc(x) ⋅ Efunc(x)
      B² = Bfunc(x) ⋅ Bfunc(x)
      vE² = E² / B²
   end

   ∇vE² = SVector{3}(ForwardDiff.gradient(E2B, X))

   bfunc(x) = normalize(Bfunc(x, t))
   ∇b̂ = ForwardDiff.jacobian(bfunc, X)
   # ∇ × b̂
   curlb = SVector{3}(∇b̂[3,2] - ∇b̂[2,3], ∇b̂[1,3] - ∇b̂[3,1], ∇b̂[2,1] - ∇b̂[1,2])
   # effective EM fields
   Eᵉ = @. E - (μ*∇B - 0.5*m*∇vE²) / q
   Bᵉ = @. B + y[4] / q2m * curlb 
   Bparᵉ = b̂ ⋅ Bᵉ # effective B field parallel to B

   dy[1] = (y[4] * Bᵉ[1] + Eᵉ[2]*b̂[3] - Eᵉ[3]*b̂[2]) / Bparᵉ
   dy[2] = (y[4] * Bᵉ[2] + Eᵉ[3]*b̂[1] - Eᵉ[1]*b̂[3]) / Bparᵉ
   dy[3] = (y[4] * Bᵉ[3] + Eᵉ[1]*b̂[2] - Eᵉ[2]*b̂[1]) / Bparᵉ
   dy[4] = q2m / Bparᵉ * Bᵉ ⋅ Eᵉ

   return
end

"1st order approximation of guiding center equations."
function trace_gc_1st!(dy, y, p::GCTuple, t)
   q, m, μ, Efunc, Bfunc = p
   q2m = q / m
   X = @view y[1:3]
   E = Efunc(X, t)
   B = Bfunc(X, t)
   b̂ = normalize(B) # unit B field at X
   u = y[4]
 
   Bmag(x) = √(Bfunc(x) ⋅ Bfunc(x))
   ∇B = SVector{3}(ForwardDiff.gradient(Bmag, X))
   Ω = q * Bmag(X) / m
 
   bfunc(x) = normalize(Bfunc(x, t))
   ∇b̂ = ForwardDiff.jacobian(bfunc, X)
   # effective EM fields
   Eᵉ = @. E - (μ*∇B) / q
   κ = ∇b̂ * b̂  # curvature
   vX = u * b̂ + b̂ × (κ*u^2 - q2m*Eᵉ) / Ω
 
   dy[1] = vX[1]
   dy[2] = vX[2]
   dy[3] = vX[3]
   dy[4] = q2m *(b̂ + u*(b̂ × κ)/Ω) ⋅ Eᵉ
 
   return
end