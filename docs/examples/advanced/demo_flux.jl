# ---
# title: Flux
# id: demo_flux
# date: 2024-04-13
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.11.4
# description: Estimation of particle flux
# ---

# This demo shows how to estimate particle flux using TestParticle.jl.
# We assume zero EM fields with constant particle velocities along the x-direction.

using TestParticle
using OrdinaryDiffEq
using StaticArrays

zeroB(x) = SA[0.0, 0.0, 0.0]
zeroE(x) = SA[0.0, 0.0, 0.0]

"Set initial conditions."
function prob_func(prob, i, repeat)
   it = (i - 1) ÷ source_flux_tp
   ir = (i - 1) % source_flux_tp
   ## initial velocity, [m/s]
   if ir < source_flux_tp ÷ 2
      v₀ = [1.0, 0.0, 0.0]
   else
      v₀ = [2.0, 0.0, 0.0]
   end
   ## initial position, [m]
   r₀ = prob.u0[1:3]
   t = (prob.tspan[1] + it, prob.tspan[2])

   prob = remake(prob; u0 = [r₀..., v₀...], tspan=t)
end

## Source flux at the origin
source_flux = 100 # [real particle / s]

param = prepare(zeroE, zeroB)
stateinit = zeros(6) # particle position and velocity to be modified
## Give particles enough time to reach steady state
tspan = (-50.0, 100.0)
const tp_per_rp = 10 # test particle per real particle
const launch_time = 20 # [s], test particle launching lasting time
const source_flux_tp = source_flux ÷ tp_per_rp # [test particles / s]
trajectories = source_flux ÷ tp_per_rp * launch_time # total number of test particles

# The number of test particles is usually smaller than the real particles. In this case, one test particle represents `tp_per_rp` real particles. Assuming steady state, the time for recording the particles is the same as the launch time. We also assume that no particles have reached the plane before we count.

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)

sols = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories)

"Estimate the particle flux through a plane x = x0 in a given time range."
function estimate_flux(sols, x0, trange=sols[1].prob.tspan)
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
flux = estimate_flux(sols, plane_loc, count_time)
println("Particle flux through plane x = $plane_loc: ", flux, " /s")

# The estimated particle flux shall match the source flux in this example. However, in general cases the particle flux through a given surface should be smaller than the source flux. We should also take the area into account.