using TestParticle, OrdinaryDiffEq, StaticArrays, Random
using TestParticle: Field, qᵢ, mᵢ, qₑ, mₑ, c, guiding_center, get_gc
using Meshes: CartesianGrid
using Test

"Initial state perturbation for EnsembleProblem."
function prob_func(prob, i, repeat)
   remake(prob, u0=rand(MersenneTwister(i))*prob.u0)
end

function prob_func_boris_mutable(prob, i, repeat)
   prob.u0[5] = i*1e5

   prob
end

function prob_func_boris_immutable(prob, i, repeat)
   prob = @views remake(prob; u0=[prob.u0[1:4]..., i*1e5, prob.u0[6]])
end

"Test boundary check method."
function isoutofdomain(u, p, t)
   if hypot(u[1], u[2], u[3]) > 0.8
      return true
   else
      return false
   end
end

@testset "TestParticle.jl" begin
   @testset "sampling" begin
      u0 = [0.0, 0.0, 0.0]
      p = 1e-9  # [Pa]
      n = 1e6  # [/cc]
      vdf = Maxwellian(u0, p, n)
      Random.seed!(1234)
      v = sample(vdf)
      @test sum(v) == 371365.50994738773
      @test startswith(repr(vdf), "Isotropic")
      B = [1.0, 0.0, 0.0] # will be normalized internally
      vdf = BiMaxwellian(B, u0, p, p, n)
      Random.seed!(1234)
      v = sample(vdf)
      @test sum(v) == -961387.4020494563
      @test startswith(repr(vdf), "BiMaxwellian")
   end

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

      grid = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position, [m]
      u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]
      tspan = (0.0, 1.0)

      param = prepare(x, y, z, E, B)
      prob = ODEProblem(trace!, stateinit, tspan, param)
      sol = solve(prob, Tsit5(); save_idxs=[1], isoutofdomain, verbose=false)
      # There are numerical differences on x86 and ARM platforms!
      @test getindex.(sol.u, 1)[end] ≈ 0.7388945226814018
      # Because the field is uniform, the order of interpolation does not matter.
      param = prepare(x, y, z, E, B; order=2)
      prob = remake(prob; p=param)
      sol = solve(prob, Tsit5(); save_idxs=[1], isoutofdomain, verbose=false)
      @test getindex.(sol.u, 1)[end] ≈ 0.7388945226814018

      param = prepare(x, y, z, E, B; order=3)
      prob = remake(prob; p=param)
      sol = solve(prob, Tsit5(); save_idxs=[1], isoutofdomain, verbose=false)
      @test getindex.(sol.u, 1)[end] ≈ 0.7388945226814018

      param = prepare(grid, E, B)
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
      tspan = (0.0, 1.0)

      param = prepare(grid, E, B)
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
      tspan = (0.0, 1.0)

      param = prepare(TestParticle.getE_dipole, TestParticle.getB_dipole)
      prob = ODEProblem(trace!, stateinit, tspan, param)
      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test guiding_center([stateinit..., 0.0], param)[1] == 1.59275e7
      @test get_gc(param) isa Function
      @test x[300] ≈ 1.2563192407332942e7 rtol=1e-6

      # static array version (results not identical with above: maybe some bugs?)
      stateinit = SA[r₀..., v₀...]

      prob = ODEProblem(trace, stateinit, tspan, param)
      sol = solve(prob, Tsit5(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test x[306] ≈ 1.2440619301099773e7 rtol=1e-5
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

      grid = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position, [m]
      u0 = [0.0, 1.0, 0.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]
      tspan = (0.0, 1.0)

      param = prepare(grid, E_field, B, F; species=Electron)
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

      grid = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
         (x[1], y[1], z[1]),
         (Δx, Δy, Δz))

      x0 = [0.0, 0.0, 0.0] # initial position, [m]
      u0 = [0.0, 1.0, 1.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]
      tspan = (0.0, 1.0)

      param = prepare(grid, E_field, B_field, F; species=Electron)
      prob = ODEProblem(trace!, stateinit, tspan, param)
      sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

      x = getindex.(sol.u, 1)
      y = getindex.(sol.u, 2)
      z = getindex.(sol.u, 3)

      @test x[end] ≈ -1.2828663442681638 && y[end] ≈ 1.5780464321537067 && z[end] ≈ 1.0
      @test guiding_center([stateinit..., 0.0], param)[1] == -0.5685630064930044

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
      # Test unsupported function form
      F_field(r, v, t) = SA[r, v, t]

      @test typeof(Field(F_field)).parameters[1] == false

      x0 = [10.0, 10.0, 0.0] # initial position, [m]
      u0 = [1e10, 0.0, 0.0] # initial velocity, [m/s]
      stateinit = [x0..., u0...]
      tspan = (0.0, 2e-7)

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
      function calc_energy(sol)
         v = hypot(sol.u[end][4:6]...)
         γ = 1/sqrt(1-(v/c)^2)
         return -(γ-1)*mₑ*c^2/qₑ
      end

      x0 = [10.0, 10.0, 0.0] # initial position, [m]
      u0 = [0.0, 0.0, 0.0] # initial velocity, [m/s]
      tspan = (0.0, 2e-7)
      stateinit = [x0..., u0...]

      param = prepare(E_field, B_field; species=Electron)
      prob = ODEProblem(trace_relativistic!, stateinit, tspan, param)
      sol = solve(prob, Vern6(); dtmax=1e-10)

      x = sol.u[end][1:3]
      # Test whether the kinetic energy [eV] of the electron
      # is equal to the electric potential energy gained.
      @test calc_energy(sol)/(x[1]-x0[1]+x[2]-x0[2]) ≈ 1e5

      prob = ODEProblem(trace_relativistic, SA[stateinit...], tspan, param)
      sol = solve(prob, Vern6(); dtmax=1e-10)
      x = sol.u[end][1:3]

      @test calc_energy(sol)/(x[1]-x0[1]+x[2]-x0[2]) ≈ 1e5
      # Tracing relativistic particle in dimensionless units
      param = prepare(xu -> SA[0.0, 0.0, 0.0], xu -> SA[0.0, 0.0, 1.0]; species=User)
      tspan = (0.0, 1.0) # 1/2π period
      stateinit = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0]
      prob = ODEProblem(trace_relativistic_normalized!, stateinit, tspan, param)
      sol = solve(prob, Vern6())
      @test sol.u[end][1] ≈ 0.46557792820784516
   end

   @testset "normalized fields" begin
      B₀ = 10e-9  # [T]
      Ω = abs(qᵢ) * B₀ / mᵢ
      t₀ = 1 / Ω  # [s]
      U₀ = 1.0    # [m/s]
      l₀ = U₀ * t₀ # [m]
      E₀ = U₀*B₀ # [V/m]

      # 3D
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      z = range(-10, 10, length=25)
      B = fill(0.0, 3, length(x), length(y), length(z)) # [B₀]
      E = fill(0.0, 3, length(x), length(y), length(z)) # [E₀]
      # This is already the normalized field; the original field is B.*B₀
      B[3,:,:,:] .= 1.0
      # E shall be normalized by E₀
      E ./= E₀

      param = prepare(x, y, z, E, B; species=User)

      x0 = [0.0, 0.0, 0.0] # initial position [l₀]
      u0 = [1.0, 0.0, 0.0] # initial velocity [v₀]
      stateinit = [x0..., u0...]
      tspan = (0.0, 0.5π) # 1/4 gyroperiod

      prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
      sol = solve(prob, Vern9(); save_idxs=[1])

      x = getindex.(sol.u, 1)

      @test length(x) == 7 && x[end] ≈ 1.0

      # 2D
      x = range(-10, 10, length=15)
      y = range(-10, 10, length=20)
      B = fill(0.0, 3, length(x), length(y)) # [B₀]
      E = fill(0.0, 3, length(x), length(y)) # [E₀]
      # This is already the normalized field; the original field is B.*B₀
      B[3,:,:] .= 1.0
      # E shall be normalized by E₀
      E ./= E₀

      Δx = x[2] - x[1]
      Δy = y[2] - y[1]

      grid = CartesianGrid((length(x)-1, length(y)-1), (x[1], y[1]), (Δx, Δy))

      x0 = [0.0, 0.0, 0.0] # initial position [l₀]
      u0 = [1.0, 0.0, 0.0] # initial velocity [v₀]
      stateinit = [x0..., u0...]
      tspan = (0.0, 0.5π) # 1/4 gyroperiod
      # periodic BC
      param = prepare(grid, E, B; species=Proton, bc=2)
      @test param[3] isa TestParticle.Field
      param = prepare(x, y, E, B; species=Proton, bc=2)
      prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
      sol = solve(prob, Tsit5(); save_idxs=[1])

      xs = getindex.(sol.u, 1)
      @test length(xs) == 9 && xs[end] ≈ 0.9999998697180689

      # Because the field is uniform, the order of interpolation does not matter.
      param = prepare(grid, E, B; order=2)
      prob = remake(prob; p=param)
      sol = solve(prob, Tsit5(); save_idxs=[1])
      xs = getindex.(sol.u, 1)
      @test length(xs) == 9 && xs[end] ≈ 0.9999998697180689

      # Because the field is uniform, the order of interpolation does not matter.
      param = prepare(grid, E, B; order=3)
      prob = remake(prob; p=param)
      sol = solve(prob, Tsit5(); save_idxs=[1])
      xs = getindex.(sol.u, 1)
      @test length(xs) == 9 && xs[end] ≈ 0.9999998697180689

      # 1D
      x = range(-10, 10, length=15)
      B = fill(0.0, 3, length(x)) # [B₀]
      E = fill(0.0, 3, length(x)) # [E₀]
      # This is already the normalized field; the original field is B.*B₀
      B[3,:] .= 1.0
      # E shall be normalized by E₀
      E ./= E₀

      x0 = [0.0, 0.0, 0.0] # initial position [l₀]
      u0 = [1.0, 0.0, 0.0] # initial velocity [v₀]
      stateinit = [x0..., u0...]
      tspan = (0.0, 0.5π) # 1/4 gyroperiod
      # periodic BC
      param = prepare(x, E, B; species=Proton, bc=3)
      prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
      sol = solve(prob, Tsit5(); save_idxs=[1])

      xs = getindex.(sol.u, 1)
      @test length(xs) == 9 && xs[end] ≈ 0.9999998697180689
   end

   @testset "Boris pusher" begin
      uniform_B(x) = SA[0.0, 0.0, 0.01]
      uniform_E(x) = SA[0.0, 0.0, 0.0]

      x0 = [0.0, 0.0, 0.0]
      v0 = [0.0, 1e5, 0.0]
      stateinit = [x0..., v0...]
      tspan = (0.0, 3e-8)
      dt = 3e-11
      param = prepare(uniform_E, uniform_B, species=Electron)
      prob = TraceProblem(stateinit, tspan, param)

      sol = TestParticle.solve(prob; dt, savestepinterval=10)[1]

      @test sol.u[end] == [-0.00010199139098074829, 3.4634030517007745e-5, 0.0,
         -62964.170425493256, -77688.56571355555, 0.0]

      t = tspan[2] / 2
      @test sol(t) == [-3.8587891411024776e-5, 5.3855910044312875e-5, 0.0,
         -93808.49725349642, 34640.52313462885, 0.0]

      prob = TraceProblem(stateinit, tspan, param; prob_func=prob_func_boris_mutable)
      trajectories = 4
      savestepinterval = 1000
      sols = TestParticle.solve(prob, EnsembleThreads(); dt, savestepinterval, trajectories)
      @test sum(s -> sum(s.u[end]), sols) ≈ -1.4065273620640622e6

      prob = TraceProblem(stateinit, tspan, param; prob_func=prob_func_boris_immutable)
      trajectories = 2
      savestepinterval = 1000
      sols = TestParticle.solve(prob; dt, savestepinterval, trajectories)
      @test sum(s -> sum(s.u[end]), sols) ≈ -421958.20861921855
   end
end

#if "makie" in ARGS
#   include("test_Makie.jl")
#end
