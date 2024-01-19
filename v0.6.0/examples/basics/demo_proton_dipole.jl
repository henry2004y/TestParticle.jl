import DisplayAs #hide
using TestParticle
using TestParticle: getB_dipole, getE_dipole, sph2cart, fieldline, mᵢ, qᵢ, c, Rₑ
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png")

# initial particle energy
Ek = 5e7 # [eV]

# initial velocity, [m/s]
v₀ = sph2cart(c*sqrt(1-1/(1+Ek*qᵢ/(mᵢ*c^2))^2), 0.0, π/4)
# initial position, [m]
r₀ = sph2cart(2.5*Rₑ, 0.0, π/2)
stateinit = [r₀..., v₀...]
# obtain field
param = prepare(getE_dipole, getB_dipole)
tspan = (0.0, 10.0)

prob = ODEProblem(trace!, stateinit, tspan, param)

sol = solve(prob, Vern9())

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
l = lines!(ax, sol)
scale!(l, invRE, invRE, invRE)

for ϕ in range(0, stop=2*π, length=10)
   lines!(fieldline(ϕ)..., color=:tomato, alpha=0.3)
end

f = DisplayAs.PNG(f) #hide

function get_energy_ratio(sol)
   vx = getindex.(sol.u, 4)
   vy = getindex.(sol.u, 5)
   vz = getindex.(sol.u, 6)

   Einit = vx[1]^2 + vy[1]^2 + vz[1]^2
   Eend = vx[end]^2 + vy[end]^2 + vz[end]^2

   (Eend - Einit) / Einit
end

function get_energy_ratio(traj::Matrix)
   Einit = traj[4,1]^2 + traj[5,1]^2 + traj[6,1]^2
   Eend = traj[4,end]^2 + traj[5,end]^2 + traj[6,end]^2
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
sol = solve(prob, Tsit5(); reltol=1e-4)

using DiffEqCallbacks

# p = (charge_mass_ratio, E, B)
dtFE(u, p, t) = 1 / (2π * abs(p[1]) * hypot(p[3](u, t)...))
cb = StepsizeLimiter(dtFE; safety_factor=1 // 10, max_step=true)

sol = solve(prob, Vern9(); callback=cb, dt=0.1) # dt=0.1 is a dummy value
get_energy_ratio(sol)

dt = 1e-4
paramBoris = BorisMethod(param)
prob = TraceProblem(stateinit, tspan, dt, paramBoris)
traj = trace_trajectory(prob)
get_energy_ratio(traj)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
