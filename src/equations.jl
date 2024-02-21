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

const FTLError = """
Particle faster than the speed of light!

If the initial velocity is slower than light and
adaptive timestepping of the solver is turned on, it
is better to set a small initial stepsize (dt) or
maximum dt for adaptive timestepping (dtmax).

More details about the keywords of initial stepsize
can be found in this documentation page:
https://diffeq.sciml.ai/stable/basics/common_solver_opts/#Stepsize-Control
"""

"""
    trace_relativistic!(dy, y, p::TPTuple, t)

ODE equations for relativistic charged particle moving in static EM field with in-place
form.
"""
function trace_relativistic!(dy, y, p::TPTuple, t)
   q2m, E, B = p

   u2 = y[4]^2 + y[5]^2 + y[6]^2
   c2 = c^2
   if u2 ≥ c2
      throw(DomainError(u2, FTLError))
   end

   γInv3 = √(1.0 - u2/c2)^3
   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   dy[1], dy[2], dy[3] = vx, vy, vz
   dy[4] = q2m*γInv3*(vy*Bz - vz*By + Ex)
   dy[5] = q2m*γInv3*(vz*Bx - vx*Bz + Ey)
   dy[6] = q2m*γInv3*(vx*By - vy*Bx + Ez)

   return
end

"""
    trace_relativistic(y, p::TPTuple, t) -> SVector{6, Float64}

ODE equations for relativistic charged particle moving in static EM field with out-of-place
form.
"""
function trace_relativistic(y, p::TPTuple, t)
   q2m, E, B = p

   u2 = y[4]^2 + y[5]^2 + y[6]^2
   c2 = c^2
   if u2 ≥ c2
      throw(DomainError(u2, FTLError))
   end

   γInv3 = √(1.0 - u2/c2)^3
   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   dx, dy, dz = vx, vy, vz
   dux = q2m*γInv3*(vy*Bz - vz*By + Ex)
   duy = q2m*γInv3*(vz*Bx - vx*Bz + Ey)
   duz = q2m*γInv3*(vx*By - vy*Bx + Ez)

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

Normalized ODE equations for relativistic charged particle moving in static EM field with
in-place form.
"""
function trace_relativistic_normalized!(dy, y, p::TPNormalizedTuple, t)
   Ω, E, B = p

   u2 = y[4]^2 + y[5]^2 + y[6]^2
   if u2 ≥ 1
      throw(DomainError(u2, FTLError))
   end

   γInv3 = √(1.0 - u2)^3
   vx, vy, vz = @view y[4:6]
   Ex, Ey, Ez = E(y, t)
   Bx, By, Bz = B(y, t)

   dy[1], dy[2], dy[3] = vx, vy, vz
   dy[4] = Ω*γInv3*(vy*Bz - vz*By + Ex)
   dy[5] = Ω*γInv3*(vz*Bx - vx*Bz + Ey)
   dy[6] = Ω*γInv3*(vx*By - vy*Bx + Ez)

   return
end