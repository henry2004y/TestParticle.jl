using TestParticle, Meshes
using Test

"Convert from spherical coordinates vector to Cartesian vector."
function sph2cart(r, ϕ, θ)
   r*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
end

"Analytic electirc field function for testing."
function getE(xu)
   [0.0, 0.0, 0.0]
end

function getB(xu)
   BMoment = [0.0, 0.0, 7.94e22]
   TestParticle.dipole(xu[1:3], BMoment)
end

@testset "TestParticle.jl" begin
   @testset "numerical field" begin
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      z = range(-10, 10, length=25)
      B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
      E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

      B[3,:,:,:] .= 10e-9
      E[3,:,:,:] .= 5e-10

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]
      Δz = z[2] - z[1]

      mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position, [m]
      u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]

      param = prepare(mesh, E, B)
      tspan = (0.0,1.0)
      sol = simulate(param, stateinit, tspan)

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test length(x) == 8 && x[end] ≈ 0.8540967226885379
   end

   @testset "analytical field" begin
      Ek = 5e7 # [eV]

      m = TestParticle.mᵢ
      q = TestParticle.qᵢ
      c = TestParticle.c
      Rₑ = TestParticle.Rₑ     

      # initial velocity, [m/s]
      v₀ = sph2cart(c*sqrt(1-1/(1+Ek*q/(m*c^2))^2), 0.0, π/4)
      # initial position, [m]
      r₀ = sph2cart(2.5*Rₑ, 0.0, π/2)
      stateinit = [r₀..., v₀...]

      param = prepare(getE, getB)
      tspan = (0.0,1.0)
      sol = simulate(param, stateinit, tspan, isAnalytic=true)

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test length(x) == 300 && x[end] ≈ 1.2592311352654776e7
   end
end
