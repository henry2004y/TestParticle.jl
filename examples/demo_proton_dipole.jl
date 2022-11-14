# Tracing charged particle in a static analytic dipole magnetic field.
# This example shows how to trace protons of a certain energy in a analytic
# Earth-like magnetic dipole field. Pyplot is used for interactively changing
# views in 3D.
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using TestParticle: getB_dipole, getE_dipole, sph2cart, Rₑ
using OrdinaryDiffEq
using PyPlot


Ek = 5e7 # [eV]

m = TestParticle.mᵢ
q = TestParticle.qᵢ
c = TestParticle.c
Rₑ = TestParticle.Rₑ     

# initial velocity, [m/s]
v₀ = sph2cart(c*sqrt(1-1/(1+Ek*q/(m*c^2))^2), 0.0, π/4)
# initial position, [m]
r₀ = sph2cart(2.5*Rₑ, 0.0, π/2)
stateinit = [r₀..., v₀...]
# obtain field
param = prepare(getE_dipole, getB_dipole)
tspan = (0.0, 2.0)

prob = ODEProblem(trace!, stateinit, tspan, param)

sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

x = getindex.(sol.u,1) / Rₑ
y = getindex.(sol.u,2) / Rₑ
z = getindex.(sol.u,3) / Rₑ

## Visualization

using3D()
fig = plt.figure()
ax = fig.gca(projection="3d")

plot(x,y,z)

for ϕ in range(0, stop=2*π, length=10)
   plot(fieldline(ϕ)..., color="r", alpha=0.3)
end