
using TestParticle
using TestParticle: getB_mirror
using OrdinaryDiffEq
using StaticArrays
using TestParticleMakie
using GLMakie
using FieldTracer
using Vlasiator: kB, rotateTensorToVectorZ
using Random
using LinearAlgebra: Symmetric
using Statistics: mean
using LinearAlgebra: ×
using RecursiveArrayTools
#using ForwardDiff: Dual

function location!(dx, x, v, p::TestParticle.TPTuple, t)
   dx .= v
end

function lorentz!(dv, x, v, p::TestParticle.TPTuple, t)
   q, m, E, B = p
   dv .= q/m*(E(x, t) + v × (B(x, t)))
end

## Obtain field

# Magnetic bottle parameters in SI units
const I1 = 1. # current in the solenoid
const N1 = 1  # number of windings
const distance = 10. # distance between solenoids
const a = 2.5 # radius of each coil
# Unfortunately there is no analytical differentiation to that special function!
function getB(xu)
   SVector{3}(getB_mirror(xu[1], xu[2], xu[3], distance, a, I1*N1))
   #SVector{3}(0.0, 0.0, 5.88e-8)
end

function getE(xu)
   SA[0.0, 0.0, 0.0]
end

function isoutofdomain(u, p, t)
   if abs(u[3]) > 6.0
      return true
   else
      return false
   end
end

struct Trajectory{T<:AbstractFloat}
   x::Vector{T}
   y::Vector{T}
   z::Vector{T}
   vx::Vector{T}
   vy::Vector{T}
   vz::Vector{T}
end

"Uniform sampling on a sphere."
function sample_unit_sphere(n::Int=1000, method::Symbol=:standard)
   x = zeros(n)
   y = zeros(n)
   z = zeros(n)

   if method == :spherical
      for i in 1:n
         theta = 2π * rand()
         phi = acos(1 - 2 * rand())
         x[i] = sin(phi) * cos(theta)
         y[i] = sin(phi) * sin(theta)
         z[i] = cos(phi)
      end
   elseif method == :standard
      i = 0
      while i < n
         x1, x2, x3 = rand(3) .* 2 .- 1
         λ = hypot(x1, x2, x3)
         if λ ≤ 1
            x[i+1] = x1 / λ
            y[i+1] = x2 / λ
            z[i+1] = x3 / λ
            i += 1
         end
      end
   end

   x, y, z
end

function sample_equator(n::Int, rmin::Float64=0.1, rmax::Float64=1.0)
   x = zeros(n)
   y = zeros(n)
   i = 0
   while i < n
      x1, x2 = (rand(2) .* 2 .- 1) .* rmax
      λ = hypot(x1, x2)
      if rmin ≤ λ ≤ rmax
         x[i+1] = x1
         y[i+1] = x2
         i += 1
      end
   end
   z = zeros(n)
   x, y, z
end

function sample_particles(n::Int=10, energy::Float64=1e-7,
   rmin::Float64=0.0, rmax::Float64=0.1)
   x, y, z = sample_equator(n, rmin, rmax)
   #x, y, z = 0.0, 1.0, 0.0
   vx, vy, vz = sample_unit_sphere(n)
   T = energy * TestParticle.qᵢ / kB
   vth = √(2kB*T/TestParticle.mᵢ)
   vx .*= vth
   vy .*= vth
   vz .*= vth
   #vx .= 0.0
   #vy .= 0.0
   #vz .= vth
   x, y, z, vx, vy, vz
end

function trace_particle(; species::TestParticle.Species=Proton, tspan=(0.0, 1.0),
   trajectories::Int=1, energy=1e-7, saveat::Float64=0.0)

   #stateinit = SA[0.0, 0.0, 5.0, 0.0, 0.0, 0.0]
   stateinit = [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]

   # Estimated max time step (better pick a location with maximum B strength!)
   Bestimate = getB(stateinit)
   Ω = TestParticle.qᵢ * hypot(Bestimate...) / TestParticle.mᵢ
   T = 2π / Ω

   param = prepare(getE, getB; species)

   #prob = ODEProblem(TestParticle.trace, stateinit, tspan, param)
   prob = DynamicalODEProblem(location!, lorentz!, stateinit[1:3], stateinit[4:6], tspan, param)

   sols = Vector{Trajectory}(undef, trajectories)

   particles = sample_particles(trajectories, energy)

   for i in 1:trajectories
      #u0 = SVector{6}([p[i] for p in particles])
      #prob = remake(prob; u0)

      u0 = [p[i] for p in particles]
      prob = remake(prob; u0=ArrayPartition(u0[1:3], u0[4:6]))

      if saveat == 0.0
         #sol = solve(prob, Tsit5(); isoutofdomain, verbose=false)
         sol = solve(prob, Tsit5(); dtmax=T/15, isoutofdomain, verbose=false)
         #sol = solve(prob, KahanLi6(); dt=0.001, isoutofdomain, verbose=true)
         #sol = solve(prob, Trapezoid(autodiff=false); dtmax=T/15, isoutofdomain, verbose=false) # bad
         #sol = solve(prob, ImplicitMidpoint(autodiff=false); dt=0.001, isoutofdomain, verbose=true) # perfect at dt=1e-3
         #sol = solve(prob, RosenbrockW6S4OS(autodiff=false); dt=0.001, isoutofdomain, verbose=true) # almost perfect at dt=1e-3
         #sol = solve(prob, SERK2(); isoutofdomain, verbose=false) # very good
         #sol = solve(prob, ESERK5(); isoutofdomain, verbose=true) # ok
         # error for ODEProblems?
         #sol = solve(prob, SSPSDIRK2(autodiff=false); dt=0.001, isoutofdomain, verbose=true) # too constrained time steps
         #sol = solve(prob, Vern6(); dtmax=T/15, isoutofdomain, verbose=false) # good
         i == 1 && println("Tmax = ", maximum(diff(sol.t)) / T)
      else
         sol = solve(prob, Tsit5();
            save_idxs=[1,2,3,4,5,6], isoutofdomain, saveat, verbose=false)
      end
      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)
      vx = getindex.(sol.u, 4)
      vy = getindex.(sol.u, 5)
      vz = getindex.(sol.u, 6)
      sols[i] = Trajectory(x, y, z, vx, vy, vz)
   end

   sols
end

function getmoments(x, y, z, vx, vy, vz, region; species::TestParticle.Species=Proton)
   nP = length(x)
   n = 0
   v = zeros(eltype(vx), 3)
   p = zeros(eltype(vx), 3, 3)
   for i in 1:nP
      if region[1] ≤ x[i] ≤ region[2] &&
         region[3] ≤ y[i] ≤ region[4] &&
         region[5] ≤ z[i] ≤ region[6]
         n += 1
         v[1] += vx[i]
         v[2] += vy[i]
         v[3] += vz[i]
         p[1,1] += vx[i]^2
         p[2,2] += vy[i]^2
         p[3,3] += vz[i]^2
         # only construct the upper matrix due to symmetry
         p[1,2] += vx[i]*vy[i]
         p[2,3] += vy[i]*vz[i]
         p[1,3] += vz[i]*vx[i]
      end
   end

   v ./= n
   if species == Proton
      p .*= TestParticle.mᵢ
   elseif species == Electron
      p .*= TestParticle.mₑ
   end

   p = Symmetric(p)

   n, v, p
end

function get_anisotropy(B::AbstractVector, P::AbstractMatrix)
   PR = rotateTensorToVectorZ(P, B)
   Paniso = 0.5f0*(PR[1,1] + PR[2,2]) / PR[3,3]
end

#####

xrange = range(-2, 2, length=20)
yrange = range(-2, 2, length=20)
zrange = range(-5, 5, length=20)

Bx, By, Bz = let x=xrange, y=yrange, z=zrange
   Bx = zeros(length(x), length(y), length(z))
   By = zeros(length(x), length(y), length(z))
   Bz = zeros(length(x), length(y), length(z))
   for k in eachindex(z), j in eachindex(y), i in eachindex(x)
      Bx[i,j,k], By[i,j,k], Bz[i,j,k] = getB([x[i], y[j], z[k]])
   end

   Bx, By, Bz
end

f = Figure()
ax = Axis3(f[1, 1];
   aspect = :data
   )

for i in 0:8
   if i == 0
      xs, ys, zs = 0.0, 0.0, 0.0
   else
      xs, ys, zs = 1.5*cos(2π*(i-1)/8), 1.5*sin(2π*(i-1)/8), 0.0
   end
   x1, y1, z1 = FieldTracer.trace(Bx, By, Bz, xs, ys, zs, xrange, yrange, zrange;
      ds=0.1, maxstep=1000)
   lines!(ax, x1, y1, z1, color=:black)
end
f

tspan = (0.0, 1000.)
energy = 1e-7
trajectories = 1

Random.seed!(1234)

sols = trace_particle(; species=Proton, tspan, trajectories, energy)

n1, v1, P1, KE1 = let
   x = [s.x[1] for s in sols]
   y = [s.y[1] for s in sols]
   z = [s.z[1] for s in sols]
   vx = [s.vx[1] for s in sols]
   vy = [s.vy[1] for s in sols]
   vz = [s.vz[1] for s in sols]
   KE1 = mean(@. √(vx^2 + vy^2 + vz^2))
   region = [-1.0, 1.0, -1.0, 1.0, -0.1, 0.1]
   n, v, P = getmoments(x, y, z, vx, vy, vz, region; species=Proton)
   n, v, P, KE1
end

n2, v2, P2, KE2 = let
   x = [s.x[end] for s in sols]
   y = [s.y[end] for s in sols]
   z = [s.z[end] for s in sols]
   vx = [s.vx[end] for s in sols]
   vy = [s.vy[end] for s in sols]
   vz = [s.vz[end] for s in sols]
   KE2 = mean(@. √(vx^2 + vy^2 + vz^2))
   region = [-2.0, 2.0, -2.0, 2.0, -6.0, 6.0]
   n, v, P = getmoments(x, y, z, vx, vy, vz, region; species=Proton)
   n, v, P, KE2
end

A1 = get_anisotropy([0.0, 0.0, 1.0], P1)
A2 = get_anisotropy([0.0, 0.0, 1.0], P2)

println("Number of particles: ", n1, " to ", n2)
# momentum conservation
println("Avg vx: ", v1[1], " to ", v2[1])
println("Avg vy: ", v1[2], " to ", v2[2])
println("Avg vz: ", v1[3], " to ", v2[3])
# energy conservation
println("Avg energy: ", KE1, " to ", KE2)
println("Energy change ratio: ", abs(KE2-KE1)/KE1)
println("Anisotropy: ", A1, " to ", A2)

for s in sols
   lines!(ax, s.x, s.y, s.z)
end
f