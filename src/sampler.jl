# Particle sampling


# Wrapper functions to avoid type piracy when adding convenience constructors

"""
     Maxwellian(args...; kw...)

Construct a `Maxwellian` distribution. Forwards to `VelocityDistributionFunctions.Maxwellian`.

     Maxwellian(u0, p, n; m=mᵢ)

Construct a `Maxwellian` distribution with bulk velocity `u0`, thermal pressure `p`, and
number density `n` in SI units. The default particle is proton.
"""
function Maxwellian end

"""
     BiMaxwellian(args...; kw...)

Construct a `BiMaxwellian` distribution. Forwards to `VelocityDistributionFunctions.BiMaxwellian`.

     BiMaxwellian(B, u0, ppar, pperp, n; m=mᵢ)

Construct a `BiMaxwellian` distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, and number density `n` in
SI units. The default particle is proton.
"""
function BiMaxwellian end

"""
     Kappa(args...; kw...)

Construct a `Kappa` distribution. Forwards to `VelocityDistributionFunctions.Kappa`.

     Kappa(u0, p, n, kappa; m=mᵢ)

Construct a `Kappa` distribution with bulk velocity `u0`, thermal pressure `p`, number density
`n`, and spectral index `kappa` in SI units. The default particle is proton.
"""
function Kappa end

"""
     BiKappa(args...; kw...)

Construct a `BiKappa` distribution. Forwards to `VelocityDistributionFunctions.BiKappa`.

     BiKappa(B, u0, ppar, pperp, n, kappa; m=mᵢ)

Construct a `BiKappa` distribution with magnetic field `B`, bulk velocity `u0`, parallel
thermal pressure `ppar`, perpendicular thermal pressure `pperp`, number density `n`, and
spectral index `kappa` in SI units. The default particle is proton.
"""
function BiKappa end

get_thermal_speed(p, n, m) = √(2 * p / (n * m))
