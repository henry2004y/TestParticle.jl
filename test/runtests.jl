using TestParticle, Meshes, OrdinaryDiffEq, StaticArrays, Random
using Test

"Initial state perturbation for EnsembleProblem."
function prob_func(prob, i, repeat)
   remake(prob, u0=rand(MersenneTwister(i))*prob.u0)
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
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace!, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test length(x) == 8 && x[end] ≈ 0.8540967226885379

      trajectories = 10
      prob = ODEProblem(trace!, stateinit, tspan, param)
      ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
      sol = solve(ensemble_prob, Tsit5(), EnsembleThreads();
         trajectories=trajectories, save_idxs=[1,2,3])

      x = getindex.(sol.u[10].u, 1)
      y = getindex.(sol.u[10].u, 2)
      z = getindex.(sol.u[10].u, 3)
      
      @test x[7] ≈ 0.09615629718624641 rtol=1e-6

      stateinit = SA[x0..., u0...]

      param = prepare(mesh, E, B)
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test length(x) == 8 && x[end] ≈ 0.8540967226885379
   end

   @testset "analytical field" begin
      # in-place version
      Ek = 5e7 # [eV]

      m = TestParticle.mᵢ
      q = TestParticle.qᵢ
      c = TestParticle.c
      Rₑ = TestParticle.Rₑ     

      # initial velocity, [m/s]
      v₀ = TestParticle.sph2cart(c*sqrt(1-1/(1+Ek*q/(m*c^2))^2), 0.0, π/4)
      # initial position, [m]
      r₀ = TestParticle.sph2cart(2.5*Rₑ, 0.0, π/2)
      stateinit = [r₀..., v₀...]

      param = prepare(TestParticle.getE_dipole, TestParticle.getB_dipole)
      tspan = (0.0, 1.0)
      
      prob = ODEProblem(trace!, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)

      @test x[300] ≈ 1.2592311352654776e7 rtol=1e-6

      # static array version (results not identical with above: maybe some bugs?)
      stateinit = SA[r₀..., v₀...]

      prob = ODEProblem(trace, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)

      @test x[306] ≈ 1.2588844644203672e7 rtol=1e-6

   end

   @testset "mixed type fields" begin
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      z = range(-10, 10, length=25)
      B = fill(0.0, 3, length(x), length(y), length(z)) # [T]

      B[3,:,:,:] .= 10e-9
      E_field(r) = SA[0, 0, 5e-10]

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]
      Δz = z[2] - z[1]

      mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position, [m]
      u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]

      param = prepare(mesh, E_field, B)
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace!, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test length(x) == 8 && x[end] ≈ 0.8540967226885379
   end
end
