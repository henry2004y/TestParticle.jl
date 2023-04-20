# ---
# title: Magnetic dipole
# id: demo_dipole
# date: 2023-04-20
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Tracing charged particle in a static analytic dipole magnetic field
# ---

# This example shows how to trace protons of a certain energy in a analytic Earth-like
# magnetic dipole field.

using JSServe: Page # hide
Page(exportable=true, offline=true) # hide

using TestParticle
using TestParticle: getB_dipole, getE_dipole, sph2cart, fieldline, mᵢ, qᵢ, c, Rₑ
using OrdinaryDiffEq
using TestParticleMakie
using WGLMakie

## initial particle energy
Ek = 5e7 # [eV]

## initial velocity, [m/s]
v₀ = sph2cart(c*sqrt(1-1/(1+Ek*qᵢ/(mᵢ*c^2))^2), 0.0, π/4)
## initial position, [m]
r₀ = sph2cart(2.5*Rₑ, 0.0, π/2)
stateinit = [r₀..., v₀...]
## obtain field
param = prepare(getE_dipole, getB_dipole)
tspan = (0.0, 2.0)

prob = ODEProblem(trace!, stateinit, tspan, param)

sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

### Visualization

f = Figure()
ax = Axis3(f[1, 1],
   title = "50 MeV Proton trajectory in Earth's dipole field",
   xlabel = "x [Re]",
   ylabel = "y [Re]",
   zlabel = "z [Re]",
   aspect = :data,
)

invRE = 1 / Rₑ
x = getindex.(sol.u, 1) .* invRE
y = getindex.(sol.u, 2) .* invRE
z = getindex.(sol.u, 3) .* invRE

l = lines!(ax, x, y, z)

for ϕ in range(0, stop=2*π, length=10)
   lines!(fieldline(ϕ)..., color=:tomato, alpha=0.3)
end

f