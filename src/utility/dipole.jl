module Dipole
# Magnetic dipole field.

import ..μ₀, ..BMoment_Earth

"Analytic electric field function for testing."
function getE(xu)
   [0.0, 0.0, 0.0]
end

"Analytic magnetic field function for testing. Return in SI unit."
function getB(xu)
   BMoment = BMoment_Earth
   dipole(xu[1:3], BMoment)
end

"Calculates the magnetic field from a dipole with magnetic moment `M` at `r`."
function dipole(rIn, M)
   x, y, z = rIn
   r = sqrt(x^2 + y^2 + z^2)
   Coef = μ₀/(4*π*r^5)

   B = [3*x^2-r^2 3*x*y     3*x*z;
        3*y*x     3*y^2-r^2 3*y*z;
        3*z*x     3*z*y     3*z^2-r^2
       ] * M * Coef
end

"""
    fieldline(L, ϕ, nP)

Creates points on one field line of the magnetic field from a dipole.
In a centered dipole magnetic field model, the path along a given L shell can be
described as
r = L*cos²λ,
where r is the radial distance (in planetary radii) to a point on the line,
λ is its co-latitude, and L is the L-shell of interest.
"""
function fieldline(ϕ::Float64, L::Float64=2.5, nP::Int=100)

   xyz = [ sph2cart(L*sin(θ)^2,ϕ,θ) for θ in range(-π,stop=π,length=nP) ]
   x = Vector{Float64}(undef,length(xyz))
   y = Vector{Float64}(undef,length(xyz))
   z = Vector{Float64}(undef,length(xyz))

   for (i, pos) in enumerate(xyz)
      x[i],y[i],z[i] = [pos...]
   end

   (x,y,z)
end

end