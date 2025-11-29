# # Particle Flux
#
# This demo shows how to estimate particle flux using TestParticle.jl.
# For the flux estimation to be accurate, we need to guarantee that
# 1. The number of test particles are enough to avoid statistical errors.
# 2. The source is properly selected. For instance, we usually sample from plasma moments of a known density, velocity and pressure. ``n*U`` gives us the source flux in units of [particles / s / m²]. If our source plane does not include all the possible sources, we will underestimate the target flux. Therefore, we need to make sure that our source plane covers all the launching possibilities, e.g. a closed sphere. Under the steady state assumption, we use the source launching time ``T`` in the denominator of calculating the flux. This time cancels out since there is also a ``T`` in the numerator.
# The target flux is then estimated as
# ```math
# J = \mathrm{h} / \mathrm{N} \oint_\mathrm{source} n\mathbf{U}\mathrm{d}S / \oint_\mathrm{target} \mathrm{d}S
# ```
# where ``\mathrm{h}`` is the number of particles that pass through the target and ``\mathrm{N}`` is the total number of test particles.
#
# In magnetosphere studies, to estimate the surface flux from ion precipitation, we can use a prescribed EM field to trace test particles originating from a closed source sphere. After a sufficiently long tracing time, each particle will either impact the surface or not. The total number flux [particles / s] is then obtained by counting all impacting particles, while the flux density [particles / s / m²] is determined by counting the impacting particles within a specific area.

using TestParticle, OrdinaryDiffEqTsit5, StaticArrays

zeroB(x) = SA[0.0, 0.0, 0.0]
zeroE(x) = SA[0.0, 0.0, 0.0]

"""
Set initial conditions.
"""
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

   prob = remake(prob; u0 = [r₀..., v₀...], tspan = t)
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
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy = false)

sols = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories)

"""
Estimate the particle flux through a plane x = x0 in a given time range.
"""
function estimate_flux_plane(sols, x0, trange = sols[1].prob.tspan)
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

# In this example, we assume zero EM fields with constant particle velocities along the x-direction. The estimated particle flux shall match the source flux in this example. However, in general cases the particle flux through a given surface should be equal or smaller than the source flux. We should also take the area into account.
#
# The second case assumes a point source at the origin. Particles are constantly isotropically launched from the source. We try to estimate the particle flux through a sphere at radius r.

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
   ## initial velocity, [m/s]
   v₀ = sample_unit_velocity_spherical()
   ## initial position, [m]
   r₀ = prob.u0[1:3]
   t = (prob.tspan[1] + it, prob.tspan[2])

   prob = remake(prob; u0 = [r₀..., v₀...], tspan = t)
end

ensemble_prob = EnsembleProblem(prob; prob_func = prob_func_iso, safetycopy = false)

sols = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories)

"""
Estimate the particle flux through a sphere with radius r = r0 in a given time range.
"""
function estimate_flux_sphere(sols, r0, trange = sols[1].prob.tspan)
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
