# Magnetic field topology in fusion research.

using SpecialFunctions: erf

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
   r = SA[x, y, z]
   cl1 = CurrentLoop(a, I1, SA[0.0, 0.0, -0.5 * distance], SA[0.0, 0.0, 1.0])
   cl2 = CurrentLoop(a, I1, SA[0.0, 0.0, 0.5 * distance], SA[0.0, 0.0, 1.0])
   B1 = getB_loop(r, cl1)
   B2 = getB_loop(r, cl2)

   return B1 + B2
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
   cl3 = CurrentLoop(b, I2, SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.0])
   B3 = getB_loop(SA[x, y, z], cl3)

   return B + B3
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
   r_vec = SA[x, y, z]

   # magnetic field of the coils
   for i in 0:15
      θ = π/16 + i*π/8 # angle between the i-th coil and the x-axis

      # Coil center and normal
      # Center at (R_major * cos(θ), R_major * sin(θ), 0)
      # R_major = b + a
      R_major = b + a
      center = SA[R_major * cos(θ), R_major * sin(θ), 0.0]

      # Normal is toroidal direction (perpendicular to poloidal plane)
      # Poloidal plane is at angle θ. Normal is (-sin(θ), cos(θ), 0)
      normal = SA[-sin(θ), cos(θ), 0.0]

      cl = CurrentLoop(a, ICoils, center, normal)
      B_coil = getB_loop(r_vec, cl)

      Bx += B_coil[1]
      By += B_coil[2]
      Bz += B_coil[3]
   end

   # magnetic field of the plasma current
   σ = a/3 # parameter of the Gauss curve
   ϕ = atan(y, x)
   # distance to centre of plasma ring
   # Plasma ring radius R_p = a + b
   R_p = a + b
   distance = √(z^2 + (x - R_p*cos(ϕ))^2 + (y - R_p*sin(ϕ))^2)

   if distance > 0.0001
      I2_r_plasma = IPlasma * erf(distance/(σ*√2))

      # Plasma is a horizontal loop at z=0, radius R_p
      cl_plasma = CurrentLoop(R_p, I2_r_plasma, SA[0.0, 0.0, 0.0], SA[0.0, 0.0, 1.0])
      B_plasma = getB_loop(r_vec, cl_plasma)

      Bx += B_plasma[1]
      By += B_plasma[2]
      Bz += B_plasma[3]
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
