# Magnetic field topology in fusion research.

using Elliptic: ellipke
using SpecialFunctions: erf

"""
    getB_mirror(x, y, z, distance, a, I1) -> Vector{Float}

Get magnetic field from a magnetic mirror generated from two coils.

# Arguments
- `x,y,z::Float`: center location in [m].
- `distance::Float`: distance between solenoids in [m].
- `a::Float`: radius of each side coil in [m].
- `I1::Float`: current in the solenoid times number of windings in side coils.
"""
function getB_mirror(x, y, z, distance, a, I1)
   r = √(x^2 + y^2) # distance from z-axis

   # 1st loop
   z₁ = z + 0.5*distance
   k2 = (4*r*a / (z₁^2 + (a + r)^2) )
   K, E = ellipke(k2)
   Bz1 = μ₀*I1 / (2π*√(z₁^2+(a+r)^2)) * ((a^2-z₁^2-r^2)/(z₁^2+(r-a)^2)*E + K)
   Br1 = μ₀*z₁*I1/(2π*r*√(z₁^2+(a+r)^2))*((z₁^2+r^2+a^2)/(z₁^2+(r-a)^2)*E - K)
   Bx1 = Br1 * x / r
   By1 = Br1 * y / r

   # 2nd loop
   z₂ = z - 0.5*distance
   k2 = (4*r*a / (z₂^2 + (a + r)^2) )
   K, E = ellipke(k2)
   Bz2 = μ₀*I1 / (2π*√(z₂^2+(a+r)^2)) * ((a^2-z₂^2-r^2)/(z₂^2+(r-a)^2)*E + K)
   Br2 = μ₀*z₂*I1/(2π*r*√(z₂^2+(a+r)^2))*((z₂^2+r^2+a^2)/(z₂^2+(r-a)^2)*E - K)
   Bx2 = Br2 * x / r
   By2 = Br2 * y / r

   # total magnetic field
   if x == 0.0 && y == 0.0
      Bx = 0.0
      By = 0.0
      Bz = Bz1 + Bz2
   else
      Bx = Bx1 + Bx2
      By = By1 + By2
      Bz = Bz1 + Bz2
   end
   [Bx, By, Bz]
end

"""
    getB_bottle(x, y, z, distance, a, b, I1, I2) -> Vector{Float}

Get magnetic field from a magnetic bottle.
Reference: https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles

# Arguments
- `x,y,z::Float`: center location in [m].
- `distance::Float`: distance between solenoids in [m].
- `a::Float`: radius of each side coil in [m].
- `b::Float`: radius of central coil in [m].
- `I1::Float`: current in the solenoid times number of windings in side coils.
- `I2::Float`: current in the central solenoid times number of windings in the
central loop.
"""
function getB_bottle(x, y, z, distance, a, b, I1, I2)
   r = √(x^2 + y^2) # distance from z-axis

   B = getB_mirror(x, y, z, distance, a, I1)

   # central loop
   z₃ = z
   k = √(4*r*b / (z₃^2 + (b+r)^2) )
   K, E = ellipke(k)
   Bz3 = μ₀*I2 / (2π*√(z₃^2+(b+r)^2)) * ((b^2-z₃^2-r^2)/(z₃^2+(r-b)^2)*E + K)
   Br3 = μ₀*z₃*I2/(2π*r*√(z₃^2+(b+r)^2))*((z₃^2+r^2+b^2)/(z₃^2+(r-b)^2)*E - K)
   Bx3 = Br3 * x / r
   By3 = Br3 * y / r

   # total magnetic field
   if x == 0.0 && y == 0.0
      B[3] += Bz3
   else
      B[1] += Bx3
      B[2] += By3
      B[3] += Bz3
   end

   B
end

"""
    getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)

Get the magnetic field from a Tokamak topology consists of 16 coils.
Original: https://github.com/BoschSamuel/Simulation-of-a-Tokamak-Fusion-Reactor/blob/master/Simulation2.m
# Arguments
- `x,y,z::Float`: location in [m].
- `a::Float`: radius of each coil in [m].
- `b::Float`: radius of central region in [m].
- `ICoil::Float`: current in the coil times number of windings.
- `IPlasma::Float`: current of the plasma?
"""
function getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
   a *= 2

   Bx, By, Bz = 0.0, 0.0, 0.0

   # magnetic field of the coils
   for i = 0:15
      θ = π/16 + i*π/8 # angle between the i-th coil and the x-axis

      if abs(sin(θ)) > 0.01
         r1_ = x/(cos(θ)-sin(θ)*tan(atan(y, x) - θ)) - a - b
         r1 = √(r1_^2 + z^2)
         z1 = ((b + a)*cos(θ) + r1_*cos(θ) - x) / sin(θ)
      else
         r1_ = x - b - a
         r1 = √(r1_^2 + z^2)
         z1 = y
      end

      k = √(4r1*a / (z1^2 + (a + r1)^2))

      K, E = ellipke(k)
      # Bz1_ is the magnetic field in the coil frame
      Bz1_ = μ₀*ICoils/(2π*√(z1^2+(a+r1)^2))*((a^2-z1^2-r1^2)/(z1^2+(r1-a)^2)*E+K)
      # Br1_ is the magnetic field in the coil frame
      Br1_ = μ₀*z1*ICoils/(2π*r1*√(z1^2+(a+r1)^2))*((z1^2+r1^2+a^2)/(z1^2+(r1-a)^2)*E-K)
      # normal coordinates
      Bx1 = -sin(θ)*Bz1_ + Br1_*r1_/r1*cos(θ)
      By1 = cos(θ)*Bz1_ + sin(θ)*Br1_*r1_/r1
      Bz1 = Br1_ * z / r1

      # add the field of a single coil to the total field
      Bx += Bx1
      By += By1
      Bz += Bz1

      if abs(Bx) < 5e-12
         Bx = 0.0
      end
      if abs(By) < 5e-12
         By = 0.0
      end
      if abs(Bz) < 5e-12
         Bz = 0.0
      end
   end

   # magnetic field of the plasma current
   σ = a/3 # parameter of the Gauss curve
   ϕ = atan(y, x)
   # distance to centre of plasma ring
   distance = √( z^2 + (x - (a + b)*cos(ϕ))^2 + (y - (a + b)*sin(ϕ))^2 )
   I2_r_plasma = IPlasma * erf(distance/(σ*√2))

   r = hypot(x, y)
   k = √(4r*(a+b)/(z^2+((a+b)+r)^2))
   K, E = ellipke(k)
   Bz_plasma = μ₀*I2_r_plasma/(2π*√(z^2+((a+b)+r)^2))*(((a+b)^2-z^2-r^2)/(z^2+(r-(a+b))^2)*E+K)
   Br_plasma = μ₀*z*I2_r_plasma/(2π*r*√(z^2+(b+r)^2))*((z^2+r^2+(a+b)^2)/(z^2+(r-(a+b))^2)*E-K)
   Bx_plasma = Br_plasma*x/r
   By_plasma = Br_plasma*y/r

   if distance > 0.0001
      Bx += Bx_plasma
      By += By_plasma
      Bz += Bz_plasma
   end

   [Bx, By, Bz]
end


"""
    getB_tokamak_profile(x, y, z, q_profile, a, R₀, Bζ0)

Reconstruct the magnetic field distribution from a safe factor(q) profile.
The formulations are from the book "Tokamak 4th Edition" by John Wesson.
# Arguments
- `x,y,z::Float`: location in [m].
- `q_profile::Function`: profile of q. The variable of this function must be the normalized radius.
- `a::Float`: minor radius [m].
- `R₀::Float`: major radius [m].
- `Bζ0::Float`: toroidal magnetic field on axis [T].
"""
function getB_tokamak_profile(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, q_profile,
                              a::AbstractFloat, R₀::AbstractFloat, Bζ0::AbstractFloat)
   R = hypot(x, y)
   r = hypot(R-R₀, z)
   if r > a
      throw(OverflowError("out of vacuum vessel"))
   end
   θ = atan(z, R-R₀)
   Bζ = Bζ0*R₀/R
   Bθ = r*Bζ/R₀/q_profile(r/a)
   ζ = atan(y, x)

   Bx = -Bζ*sin(ζ) - Bθ*sin(θ)*cos(ζ)
   By = Bζ*cos(ζ) - Bθ*sin(θ)*sin(ζ)
   Bz = Bθ*cos(θ)

   [Bx, By, Bz]
end