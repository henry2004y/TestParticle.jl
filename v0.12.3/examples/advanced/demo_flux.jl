using TestParticle
using OrdinaryDiffEq
using StaticArrays

zeroB(x) = SA[0.0, 0.0, 0.0]
zeroE(x) = SA[0.0, 0.0, 0.0]

"Set initial conditions."
function prob_func(prob, i, repeat)
   it = (i - 1) ÷ source_flux_tp
   ir = (i - 1) % source_flux_tp
   # initial velocity, [m/s]
   if ir < source_flux_tp ÷ 2
      v₀ = [1.0, 0.0, 0.0]
   else
      v₀ = [2.0, 0.0, 0.0]
   end
   # initial position, [m]
   r₀ = prob.u0[1:3]
   t = (prob.tspan[1] + it, prob.tspan[2])

   prob = remake(prob; u0 = [r₀..., v₀...], tspan=t)
end

# Source flux at the origin
source_flux = 100 # [real particle / s]

param = prepare(zeroE, zeroB)
stateinit = zeros(6) # particle position and velocity to be modified
# Give particles enough time to reach steady state
tspan = (-50.0, 100.0)
const tp_per_rp = 10 # test particle per real particle
const launch_time = 20 # [s], test particle launching lasting time
const source_flux_tp = source_flux ÷ tp_per_rp # [test particles / s]
trajectories = source_flux ÷ tp_per_rp * launch_time # total number of test particles

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)

sols = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories)

"Estimate the particle flux through a plane x = x0 in a given time range."
function estimate_flux_plane(sols, x0, trange=sols[1].prob.tspan)
   count = 0
   for sol in sols
      xs = [sol(t)[1] for t in trange]
      if xs[1] ≤ x0 ≤ xs[2] || xs[1] ≥ x0 ≥ xs[2]
         count += 1
      end
   end
   count *= tp_per_rp
   count /= launch_time
end

plane_loc = 100.0 # [m]
count_time = (0.0, tspan[2]) # [s]
flux = estimate_flux_plane(sols, plane_loc, count_time)
println("Example 1:")
println("Particle flux through plane x = $plane_loc [m]: ", flux, " /s")

function sample_unit_velocity_spherical()
   ϕ = 2 * π * rand()    # Azimuthal angle in [0, 2π)
   cosθ = 2 * rand() - 1 # Uniformly sample cos(θ) in [-1, 1]
   sinθ = sqrt(1 - cosθ^2)
   x = sinθ * cos(ϕ)
   y = sinθ * sin(ϕ)
   z = cosθ
   return SA[x, y, z]
end

function prob_func_iso(prob, i, repeat)
   it = (i - 1) ÷ source_flux_tp
   # initial velocity, [m/s]
   v₀ = sample_unit_velocity_spherical()
   # initial position, [m]
   r₀ = prob.u0[1:3]
   t = (prob.tspan[1] + it, prob.tspan[2])

   prob = remake(prob; u0 = [r₀..., v₀...], tspan=t)
end

ensemble_prob = EnsembleProblem(prob; prob_func=prob_func_iso, safetycopy=false)

sols = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories)

"Estimate the particle flux through a sphere with radius r = r0 in a given time range."
function estimate_flux_sphere(sols, r0, trange=sols[1].prob.tspan)
   count = 0
   for sol in sols
      rs = [hypot(sol(t)[1:3]...) for t in trange]
      if rs[1] ≤ r0 ≤ rs[2] || rs[1] ≥ r0 ≥ rs[2]
         count += 1
      end
   end
   count *= tp_per_rp
   count /= launch_time
end

r = 100.0 # [m]
area = 4π*r^2 # [m²]
count_time = (0.0, tspan[2]) # [s]
flux = estimate_flux_sphere(sols, r, count_time)
println("Example 2:")
println("Particle flux through sphere r = $r [m]: ", flux, " /s")
println("Particle flux density through sphere r = $r [m]: ", flux / area, " /(s * m²)")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
