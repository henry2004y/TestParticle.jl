# Magnetic field topology in fusion research.

using Elliptic: ellipke
using SpecialFunctions: erf

#TODO Generalize the geometry
struct Currentloop{T <: AbstractFloat, Vl, Vn}
   "radius of coil loop [m]"
   a::T
   "steady current [A]"
   I::T
   "location of loop center [m]"
   location::Vl
   "unit orientation of the loop (currently only support ẑ)"
   normal::Vn
end

"""
     getB_current_loop(x, y, z, cl::Currentloop) -> StaticVector{Float64, 3}

Get magnetic field at `[x, y, z]` from a magnetic mirror generated from two coils.

# Arguments

  - `x,y,z::Float`: particle coordinates in [m].
  - `distance::Float`: distance between solenoids in [m].
  - `a::Float`: radius of each side coil in [m].
  - `I1::Float`: current in the solenoid times number of windings in side coils.
"""
function getB_current_loop(x, y, z, cl::Currentloop)
   (; a, I, location) = cl
   r = √((x - location[1])^2 + (y - location[2])^2) # distance from z-axis
   z₁ = z - location[3]
   K, E = ellipke(4*r*a / (z₁^2 + (a + r)^2))
   Bz = μ₀*I / (2π*√(z₁^2 + (a+r)^2)) * ((a^2 - z₁^2 - r^2)/(z₁^2 + (r-a)^2)*E + K)
   Br = μ₀*z₁*I / (2π*r*√(z₁^2 + (a+r)^2))*((z₁^2 + r^2 + a^2)/(z₁^2 + (r-a)^2)*E - K)
   Bx = Br * x / r
   By = Br * y / r

   SA[Bx, By, Bz]
end

"""
     getB_mirror(x, y, z, distance, a, I1) -> StaticVector{Float64, 3}

Get magnetic field at `[x, y, z]` from a magnetic mirror generated from two coils.

# Arguments

  - `x,y,z::Float`: particle coordinates in [m].
  - `distance::Float`: distance between solenoids in [m].
  - `a::Float`: radius of each side coil in [m].
  - `I1::Float`: current in the solenoid times number of windings in side coils.
"""
function getB_mirror(x, y, z, distance, a, I1)
   cl1 = Currentloop(a, I1, SA[0.0, 0.0, -0.5 * distance], SA[0.0, 0.0, 1.0])
   cl2 = Currentloop(a, I1, SA[0.0, 0.0, 0.5 * distance], SA[0.0, 0.0, 1.0])
   B1 = getB_current_loop(x, y, z, cl1)
   B2 = getB_current_loop(x, y, z, cl2)
   # total magnetic field
   if x == 0.0 && y == 0.0
      B = SA[0.0, 0.0, B1[3] + B2[3]]
   else
      B = B1 + B2
   end

   B
end

"""
     getB_bottle(x, y, z, distance, a, b, I1, I2) -> StaticVector{Float64, 3}

Get magnetic field from a magnetic bottle.
Reference: [wiki](https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles)

# Arguments

  - `x,y,z::Float`: particle coordinates in [m].
  - `distance::Float`: distance between solenoids in [m].
  - `a::Float`: radius of each side coil in [m].
  - `b::Float`: radius of central coil in [m].
  - `I1::Float`: current in the solenoid times number of windings in side coils in [A].
  - `I2::Float`: current in the central solenoid times number of windings in the
    central loop in [A].
"""
function getB_bottle(x, y, z, distance, a, b, I1, I2)
   B = getB_mirror(x, y, z, distance, a, I1)
   # Central loop
   cl3 = Currentloop(b, I2, SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.0])
   B3 = getB_current_loop(x, y, z, cl3)
   # total magnetic field
   if x == 0.0 && y == 0.0
      B = SA[0.0, 0.0, B[3] + B3[3]]
   else
      B = B + B3
   end

   B
end

"""
     getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma) -> StaticVector{Float64, 3}

Get the magnetic field from a Tokamak topology consists of 16 coils.
Original: [Tokamak-Fusion-Reactor](https://github.com/BoschSamuel/Simulation-of-a-Tokamak-Fusion-Reactor/blob/master/Simulation2.m)

# Arguments

  - `x,y,z::Float`: location in [m].
  - `a::Float`: radius of each coil in [m].
  - `b::Float`: radius of central region in [m].
  - `ICoil::Float`: current in the coil times number of windings in [A].
  - `IPlasma::Float`: current of the plasma in [A].
"""
function getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma)
   a *= 2

   Bx, By, Bz = 0.0, 0.0, 0.0

   # magnetic field of the coils
   for i in 0:15
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
   distance = √(z^2 + (x - (a + b)*cos(ϕ))^2 + (y - (a + b)*sin(ϕ))^2)
   I2_r_plasma = IPlasma * erf(distance/(σ*√2))

   r = √(x^2 + y^2)
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

   SA[Bx, By, Bz]
end

"""
     getB_tokamak_profile(x, y, z, q_profile, a, R₀, Bζ0) -> StaticVector{Float64, 3}

Reconstruct the magnetic field distribution from a safe factor(q) profile.
Reference: Tokamak, 4th Edition, John Wesson.

# Arguments

  - `x,y,z::Float`: location in [m].
  - `q_profile::Function`: profile of q. The variable of this function must be the normalized radius.
  - `a::Float`: minor radius [m].
  - `R₀::Float`: major radius [m].
  - `Bζ0::Float`: toroidal magnetic field on axis [T].
"""
function getB_tokamak_profile(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
      q_profile, a::AbstractFloat, R₀::AbstractFloat, Bζ0::AbstractFloat)
   R = √(x^2 + y^2)
   r = √((R - R₀)^2 + z^2)
   if r > a
      throw(OverflowError("out of vacuum vessel"))
   end
   θ = atan(z, R - R₀)
   Bζ = Bζ0 * R₀ / R
   Bθ = r * Bζ / R₀ / q_profile(r/a)
   ζ = atan(y, x)

   Bx = -Bζ*sin(ζ) - Bθ*sin(θ)*cos(ζ)
   By = Bζ*cos(ζ) - Bθ*sin(θ)*sin(ζ)
   Bz = Bθ*cos(θ)

   SA[Bx, By, Bz]
end
"""
     getB_zpinch(x, y, z, I, a) -> StaticVector{Float64, 3}

Get magnetic field from a Z-pinch configuration.
Reference: [Z-pinch](https://en.wikipedia.org/wiki/Z-pinch)

# Arguments

  - `x,y,z::Float`: particle coordinates in [m].
  - `I::Float`: current in the wire [A].
  - `a::Float`: radius of the wire [m].
"""
function getB_zpinch(x, y, z, I, a)
   r = hypot(x, y)
   if r < a
      factor = μ₀ * I / (2π * a^2)
      Bx = -factor * y
      By = factor * x
   else
      factor = μ₀ * I / (2π * r^2)
      Bx = -factor * y
      By = factor * x
   end

   SA[Bx, By, 0.0]
end
