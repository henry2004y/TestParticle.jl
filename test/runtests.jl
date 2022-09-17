using TestParticle, Meshes, OrdinaryDiffEq, StaticArrays, Random
using TestParticle: Field, qₑ, mₑ, c
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

      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test length(x) == 8 && x[end] ≈ 0.8540967226885379

      trajectories = 10
      prob = ODEProblem(trace!, stateinit, tspan, param)
      ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
      sol = solve(ensemble_prob, Tsit5(), EnsembleThreads();
         trajectories=trajectories, save_idxs=[1])

      x = getindex.(sol.u[10].u, 1)
      
      @test x[7] ≈ 0.09615629718624641 rtol=1e-6

      stateinit = SA[x0..., u0...]

      param = prepare(mesh, E, B)
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

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

      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test x[300] ≈ 1.2592311352654776e7 rtol=1e-6

      # static array version (results not identical with above: maybe some bugs?)
      stateinit = SA[r₀..., v₀...]

      prob = ODEProblem(trace, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test x[306] ≈ 1.2588844644203672e7 rtol=1e-6

   end

   @testset "mixed type fields" begin
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      z = range(-10, 10, length=25)
      B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
      F = fill(0.0, 3, length(x), length(y), length(z)) # [N]

      B[3,:,:,:] .= 1e-11
      E_field(r) = SA[0, 5e-11, 0]  # [V/M]
      F[1,:,:,:] .= 9.10938356e-42

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]
      Δz = z[2] - z[1]

      mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position, [m]
      u0 = [0.0, 1.0, 0.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]

      param = prepare(mesh, E_field, B, F; species=Electron)
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace!, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1,2])

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)

      @test x[end] ≈ 1.5324506965560782 && y[end] ≈ -2.8156470047903706
   end

   @testset "time-independent fields" begin
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      z = range(-10, 10, length=25)
      F = fill(0.0, 3, length(x), length(y), length(z)) # [N]

      B_field(r, t) = SA[0, 0, 1e-11*cos(2π*t)] # [T]
      E_field(r, t) = SA[5e-11*sin(2π*t), 0, 0]  # [V/M]
      F[2,:,:,:] .= 9.10938356e-42

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]
      Δz = z[2] - z[1]

      mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position, [m]
      u0 = [0.0, 1.0, 1.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]

      param = prepare(mesh, E_field, B_field, F; species=Electron)
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace!, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test x[end] ≈ -1.2828663442681638 && y[end] ≈ 1.5780464321537067 && z[end] ≈ 1.0

      F_field(r) = SA[0, 9.10938356e-42, 0] # [N]

      param = prepare(E_field, B_field, F_field; species=Electron)
      _, _, _, _, F = param

      @test F(x0)[2] ≈ 9.10938356e-42

      stateinit = SA[x0..., u0...]

      prob = ODEProblem(trace, stateinit, tspan, param)
      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test x[end] ≈ -1.2828663442681638 && y[end] ≈ 1.5780464321537067 && z[end] ≈ 1.0
   end

   @testset "Exceptions" begin
      E_field(r, t) = SA[5e-11*sin(2π*t), 0, 0] 
      E = Field(E_field)

      @test_throws ArgumentError E([0, 0, 0])

      F_field(r, v, t) = SA[r, v, t]

      @test_throws ArgumentError Field(F_field)

      x0 = [10.0, 10.0, 0.0] # initial position, [m]
      u0 = [1e10, 0.0, 0.0] # initial velocity, [m/s]
      tspan = (0.0, 2e-7)
      stateinit = [x0..., u0...]
      param = prepare(E_field, E_field; species=Electron)

      prob = ODEProblem(trace_relativistic!, stateinit, tspan, param)
      @test_throws DomainError solve(prob, Tsit5())

      prob = ODEProblem(trace_relativistic, SA[stateinit...], tspan, param)
      @test_throws DomainError solve(prob, Tsit5())
   end

   @testset "relativistic particle" begin
      B_field(xu) = SA[0, 0, 0.01]
      function E_field(xu)
         Ex = 0<=xu[1]<=100 ? -1e5 : 0.0
         Ey = 0<=xu[2]<=100 ? -1e5 : 0.0
         return SA[Ex, Ey, 0.0]
      end

      # calculate the energy [eV] of a electron
      function cal_energy(sol)
         v = sol.u[end][4:6]
         v = hypot(v...)
         γ = 1/sqrt(1-(v/c)^2)
         return -(γ-1)*mₑ*c^2/qₑ
      end

      x0 = [10.0, 10.0, 0.0] # initial position, [m]
      u0 = [0.0, 0.0, 0.0] # initial velocity, [m/s]
      tspan = (0.0, 2e-7)
      stateinit = [x0..., u0...]
      param = prepare(E_field, B_field; species=Electron)

      prob = ODEProblem(trace_relativistic!, stateinit, tspan, param)
      sol = solve(prob, Vern6(); dtmax=1e-10, save_idxs=[1,2,3,4,5,6])
      x = sol.u[end][1:3]

      # Test whether the kinetic energy [eV] of the electron
      # is equal to the electric potential energy gained.
      @test cal_energy(sol)/(x[1]-x0[1]+x[2]-x0[2]) ≈ 1e5

      prob = ODEProblem(trace_relativistic, SA[stateinit...], tspan, param)
      sol = solve(prob, Vern6(); dtmax=1e-10, save_idxs=[1,2,3,4,5,6])
      x = sol.u[end][1:3]

      @test cal_energy(sol)/(x[1]-x0[1]+x[2]-x0[2]) ≈ 1e5
   end

   @testset "normalized field" begin
      # Basic units: length l₀, time t₀, magnetic field B₀
      # Derived units: velocity v₀ = l₀/t₀, electric field E₀ = v₀B₀

      # 3D
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      z = range(-10, 10, length=25)
      B = fill(0.0, 3, length(x), length(y), length(z)) # [B₀]
      E = fill(0.0, 3, length(x), length(y), length(z)) # [E₀]

      B₀ = 10e-9

      B[3,:,:,:] .= 1.0

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]
      Δz = z[2] - z[1]

      mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position [l₀]
      u0 = [1.0, 0.0, 0.0] # initial velocity [v₀]
      stateinit = [x0..., u0...]

      param = prepare(mesh, E, B, B₀; species=Proton)
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace_normalized!, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test length(x) == 8 && x[end] ≈ 0.8540967226941674

      # 2D
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      B = fill(0.0, 3, length(x), length(y)) # [B₀]
      E = fill(0.0, 3, length(x), length(y)) # [E₀]

      B₀ = 10e-9

      B[3,:,:] .= 1.0

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]

      mesh = CartesianGrid((length(x)-1, length(y)-1), (x[1], y[1]), (Δx, Δy))

      x0 = [0.0, 0.0] # initial position [l₀]
      u0 = [1.0, 0.0, 0.0] # initial velocity [v₀]
      stateinit = [x0..., u0...]

      param = prepare(mesh, E, B, B₀; species=Proton)
      tspan = (0.0, 1.0)

      prob = ODEProblem(trace2d_normalized!, stateinit, tspan, param)

      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test length(x) == 8 && x[end] ≈ 0.8540967195469715
   end
end
