# Tracing multiple charged particles in a static EM field.
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using OrdinaryDiffEq
using GLMakie
using TestParticleMakie

struct Trajectory{T<:AbstractFloat}
   x::Vector{T}
   y::Vector{T}
   z::Vector{T}
end

function trace(x, y, z, E, B; trajectories::Int=10)
   # Initialize particles
   x0 = [0.0, 0.0, 0.0] # initial position, [m]
   u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
   stateinit = [x0..., u0...]
   
   param = prepare(x, y, z, E, B, species=Electron)
   tspan = (0.0, 15.0)

   prob = ODEProblem(trace!, stateinit, tspan, param)
   # Trajectory container
   sols = Vector{Trajectory}(undef, trajectories)

   for i in 1:trajectories
      prob = remake(prob; u0=rand()*prob.u0)

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])
      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)
      sols[i] = Trajectory(x, y, z)
   end

   return sols
end

## Initialize grid and field

x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
z = range(-10, 10, length=25)
B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

B[3,:,:,:] .= 1e-11
E[3,:,:,:] .= 5e-13

trajectories = 10

## Solve for the trajectories

sols = trace(x, y, z, E, B; trajectories)

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Particle trajectories, EGI, t=1298s",
   #xlabel = "X",
   #ylabel = "Y",
   #zlabel = "Z",
   aspect = :data,
)

for i in eachindex(sols)
   lines!(sols[i].x, sols[i].y, sols[i].z, label="$i")
end

f