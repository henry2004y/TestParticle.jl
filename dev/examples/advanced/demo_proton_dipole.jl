import DisplayAs #hide
using TestParticle
using TestParticle: getB_dipole, getE_dipole, sph2cart, dipole_fieldline, mᵢ, qᵢ, c, Rₑ
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# Initial condition
stateinit = let
   # Initial particle energy
   Ek = 5e7 # [eV]
   # initial velocity, [m/s]
   v₀ = sph2cart(c*sqrt(1-1/(1+Ek*qᵢ/(mᵢ*c^2))^2), 0.0, π/4)
   # initial position, [m]
   r₀ = sph2cart(2.5*Rₑ, 0.0, π/2)
   [r₀..., v₀...]
end
# obtain field
param = prepare(getE_dipole, getB_dipole)
tspan = (0.0, 10.0)

prob = ODEProblem(trace!, stateinit, tspan, param)

sol = solve(prob, Vern9())

### Visualization

f = Figure(fontsize=18)
##ax = Axis3(f[1, 1],
#   title = "50 MeV Proton trajectory in Earth's dipole field",
#   xlabel = "x [Re]",
#   ylabel = "y [Re]",
#   zlabel = "z [Re]",
#   aspect = :data,
#   limits = (-2.5, 2.5, -2.5, 2.5, -1, 1)
##)
ax = LScene(f[1,1])
invRE = 1 / Rₑ

l = lines!(ax, sol, idxs=(1, 2, 3))
# In Makie 0.21.11, scene scaling has no effect on Axis3.
scale!(ax.scene, invRE, invRE, invRE)

for ϕ in range(0, stop=2*π, length=10)
   lines!(dipole_fieldline(ϕ).*Rₑ..., color=:tomato, alpha=0.3)
end

f = DisplayAs.PNG(f) #hide

function get_energy_ratio(sol)
   vx = @view sol[4,:]
   vy = @view sol[5,:]
   vz = @view sol[6,:]

   Einit = vx[1]^2 + vy[1]^2 + vz[1]^2
   Eend = vx[end]^2 + vy[end]^2 + vz[end]^2

   (Eend - Einit) / Einit
end

# `ImplicitMidpoint()` requires a fixed time step.
sol = solve(prob, ImplicitMidpoint(); dt=1e-3)
get_energy_ratio(sol)

sol = solve(prob, ImplicitMidpoint(); dt=1e-4)
get_energy_ratio(sol)

sol = solve(prob, Vern9())
get_energy_ratio(sol)

sol = solve(prob, Trapezoid())
get_energy_ratio(sol)

sol = solve(prob, Vern6())
get_energy_ratio(sol)

sol = solve(prob, Tsit5())
get_energy_ratio(sol)

# This is roughly equivalent in accuracy and performance with Vern9() and `reltol=1e-3` (default)
sol = solve(prob, Tsit5(); reltol=1e-4);

using DiffEqCallbacks

# p = (charge_mass_ratio, E, B)
dtFE(u, p, t) = 2π / (abs(p[1]) * sqrt(sum(x -> x^2, p[3](u, t))))

cb = StepsizeLimiter(dtFE; safety_factor=1 // 10, max_step=true)

sol = solve(prob, Vern9(); callback=cb, dt=0.1) # dt=0.1 is a dummy value
get_energy_ratio(sol)

dt = 1e-4
prob = TraceProblem(stateinit, tspan, param)
sol = TestParticle.solve(prob; dt)[1]
get_energy_ratio(sol)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
