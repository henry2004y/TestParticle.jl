# ---
# title: Multiple particles
# id: demo_multiple
# date: 2023-04-20
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Tracing multiple charged particles in a static EM field
# ---

using JSServe: Page # hide
Page(exportable=true, offline=true) # hide

using TestParticle
using OrdinaryDiffEq
using TestParticleMakie
using WGLMakie
using Random

## For reproducible results
Random.seed!(1234)

function trace(x, y, z, E, B; trajectories::Int=10)
   ## Initialize particles
   x0 = [0.0, 0.0, 0.0] # initial position, [m]
   u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
   stateinit = [x0..., u0...]

   param = prepare(x, y, z, E, B, species=Electron)
   tspan = (0.0, 15.0)

   prob = ODEProblem(trace!, stateinit, tspan, param)

   sols = Vector{ODESolution}(undef, trajectories)
   ## Sample from a Maxwellian with bulk speed 0 and thermal speed 1.0
   vdf = Maxwellian([0.0, 0.0, 0.0], 1.0)
   v = sample(vdf, trajectories)

   for i in 1:trajectories
      prob = remake(prob; u0=[x0..., v[:,i]...])

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])
      sols[i] = sol
   end

   return sols
end

### Initialize grid and field

x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
z = range(-10, 10, length=25)
B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

B[3,:,:,:] .= 1e-11
E[3,:,:,:] .= 5e-13

trajectories = 10

### Solve for the trajectories

sols = trace(x, y, z, E, B; trajectories)

### Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Particle trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data,
)

for i in eachindex(sols)
   lines!(ax, sols[i], label="$i")
end

f