"Analytic electric field function for testing."
function getE(xu)
   [0.0, 0.0, 0.0]
end

"Analytic magnetic field function for testing. Return in SI unit."
function getB(xu)
   BMoment = [0.0, 0.0, 7.94e22] # [V*s/(A*m)]
   dipole(xu[1:3], BMoment)
end

"Calculates the magnetic field from a dipole with magnetic moment `M` at `r`."
function dipole(rIn, M)
   x, y, z = rIn
   r = sqrt(x^2 + y^2 + z^2)
   Coef = TestParticle.μ₀/(4*π*r^5)

   B = [3*x^2-r^2 3*x*y     3*x*z;
        3*y*x     3*y^2-r^2 3*y*z;
        3*z*x     3*z*y     3*z^2-r^2
       ] * M * Coef
end

"Convert from spherical to Cartesian coordinates vector."
function sph2cart(r, ϕ, θ)
   r*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end