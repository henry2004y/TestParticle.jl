# # Multiple Particles

import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner
using Random
using CairoMakie
CairoMakie.activate!(type = "png") #hide

## For reproducible results
Random.seed!(1234)

function trace(x, y, z, E, B; trajectories::Int = 10)
   ## Initialize particles
   x0 = [0.0, 0.0, 0.0] # initial position, [m]
   u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
   stateinit = [x0..., u0...]

   param = prepare(x, y, z, E, B, species = Electron)
   tspan = (0.0, 15.0)

   prob = ODEProblem(trace!, stateinit, tspan, param)

   sols = Vector{ODESolution}(undef, trajectories)
   ## Sample from a Maxwellian with bulk speed 0 and thermal speed 1.0
   vdf = Maxwellian([0.0, 0.0, 0.0], 1.0)
   v = [sample(vdf) for _ in 1:trajectories]

   for i in 1:trajectories
      #prob = remake(prob; u0=[x0..., v[:,i]...])
      prob.u0[4:6] = v[i]

      sol = solve(prob, Vern9())
      sols[i] = sol
   end

   return sols
end

### Initialize grid and field
x = range(-10, 10, length = 15)
y = range(-10, 10, length = 20)
z = range(-10, 10, length = 25)

B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]
B[3, :, :, :] .= 1e-11
E[3, :, :, :] .= 5e-13

trajectories = 4

### Solve for the trajectories

sols = trace(x, y, z, E, B; trajectories)

### Visualization
f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Particle trajectories",
   xlabel = "X [m]",
   ylabel = "Y [m]",
   zlabel = "Z [m]",
   aspect = :data
)

for i in eachindex(sols)
   lines!(ax, sols[i], idxs = (1, 2, 3), color = Makie.wong_colors()[i], label = "$i")
end

f = DisplayAs.PNG(f) #hide
