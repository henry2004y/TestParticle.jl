# # Particle Flux
#
# This demo shows how to estimate particle flux using TestParticle.jl.
# For the flux estimation to be accurate, we need to guarantee that
# 1. The number of test particles is sufficient to avoid statistical errors.
# 2. The source is properly selected. For instance, we usually sample from plasma moments of a known density, velocity and pressure. ``n \cdot U`` gives us the source flux in units of [particles / s / m²]. If our source plane does not include all the possible sources, we will underestimate the target flux. Therefore, we need to make sure that our source plane covers all the launching possibilities, e.g. a closed sphere.
#
# Under the steady state assumption, we can simulate a representative set of $N$ particles launched simultaneously. If the total physical source rate is $R$ [particles / s], each simulated particle represents a contribution of $w = R/N$ [particles / s].
#
# The target flux (rate) is then estimated as
# ```math
# J = \sum_{i=1}^{h} w = h \cdot \frac{R}{N}
# ```
# where $h$ is the number of particles that impact the target.
#
# The flux density [particles / s / m²] is obtained by dividing $J$ by the target area.
#
# In magnetosphere studies, to estimate the surface flux from ion precipitation, we can use a prescribed EM field to trace test particles originating from a closed source sphere. After a sufficiently long tracing time, each particle will either impact the surface or not.

using TestParticle, OrdinaryDiffEqTsit5, StaticArrays

zeroB(x) = SA[0.0, 0.0, 0.0]
zeroE(x) = SA[0.0, 0.0, 0.0]

"Set initial conditions."
function prob_func(prob, i, repeat)
   ## initial velocity, [m/s]
   ## 50% v=1.0, 50% v=2.0
   if i % 2 == 1
      v₀ = [1.0, 0.0, 0.0]
   else
      v₀ = [2.0, 0.0, 0.0]
   end
   ## initial position, [m]
   r₀ = prob.u0[1:3]

   prob = remake(prob; u0 = [r₀..., v₀...])
end

## Source flux at the origin
const source_flux = 100.0 # [real particle / s]
const n_particles = 1000 # number of test particles
const weight = source_flux / n_particles # [particles / s per test particle]

param = prepare(zeroE, zeroB)
stateinit = zeros(6) # particle position and velocity to be modified
tspan = (0.0, 110.0) # Give particles enough time to reach the plane (x=100 with v>=1)

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy = false)

sols = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories = n_particles)

"Estimate the particle flux through a plane x = x0."
function estimate_flux_plane(sols, x0)
   count = sum(sols.u) do sol
      ## Check if particle crossed x0
      ## Assuming monotonic motion along x, check start and end
      x_start = sol(sol.prob.tspan[1])[1]
      x_end = sol(sol.prob.tspan[2])[1]

      (x_start < x0 <= x_end) || (x_start > x0 >= x_end)
   end
   return count * weight
end

plane_loc = 100.0 # [m]
flux = estimate_flux_plane(sols, plane_loc)
println("Example 1:")
println("Particle flux through plane x = $plane_loc [m]: ", flux, " /s")

# In this example, we assume zero EM fields with constant particle velocities along the x-direction. The estimated particle flux shall match the source flux in this example.
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
   ## initial velocity, [m/s]
   v₀ = sample_unit_velocity_spherical()
   ## initial position, [m]
   r₀ = prob.u0[1:3]

   prob = remake(prob; u0 = [r₀..., v₀...])
end

ensemble_prob = EnsembleProblem(prob; prob_func = prob_func_iso, safetycopy = false)

sols = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories = n_particles)

"Estimate the particle flux through a sphere with radius r = r0."
function estimate_flux_sphere(sols, r0)
   count = sum(sols.u) do sol
      u_start = sol(sol.prob.tspan[1])
      r_start = hypot(u_start[1], u_start[2], u_start[3])
      u_end = sol(sol.prob.tspan[2])
      r_end = hypot(u_end[1], u_end[2], u_end[3])

      (r_start < r0 <= r_end) || (r_start > r0 >= r_end)
   end
   return count * weight
end

r = 100.0 # [m]
area = 4π*r^2 # [m²]
flux = estimate_flux_sphere(sols, r)
println("Example 2:")
println("Particle flux through sphere r = $r [m]: ", flux, " /s")
println("Particle flux density through sphere r = $r [m]: ", flux / area, " /(s * m²)")
