
"""

Variable `(u, x, y, z)`, where `u` is the velocity along the magnetic field.
"""
function trace_gc(y, p::TPTuple, t)
   q2m, E, B = p
   b̂ = normalize(B)
   # for simplicity
   Eᵉ = E
   Bᵉ = B
   Bparᵉ = b̂ ⋅ Bᵉ # parallel effective B field

   du = q2m / Bparᵉ * Bᵉ ⋅ Eᵉ
   dX = (y[1] * Bᵉ + Eᵉ × b̂) / Bparᵉ

   #vx, vy, vz = @view y[2:4]
   #Ex, Ey, Ez = E(y, t)
   #Bx, By, Bz = B(y, t)

   #dx, dy, dz = vx, vy, vz
   # q/m*(E + v × B)
   #dux = q2m*(vy*Bz - vz*By + Ex)
   #duy = q2m*(vz*Bx - vx*Bz + Ey)
   #duz = q2m*(vx*By - vy*Bx + Ez)
   SVector{4}(du, dX[1], dX[2], dX[3])
end