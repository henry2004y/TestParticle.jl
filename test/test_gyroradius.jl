using Test
using StaticArrays
using TestParticle
using OrdinaryDiffEq

@testset "Gyroradius Utility" begin
   @testset "Non-relativistic" begin
      # Uniform B field
      B0 = 1e-8
      B_func(x, t) = SA[0.0, 0.0, B0]
      E_func(x, t) = SA[0.0, 0.0, 0.0]

      # Particle parameters
      q = TestParticle.qᵢ
      m = TestParticle.mᵢ

      # Initial state: v perpendicular to B
      v_perp = 1e5
      x0 = SA[0.0, 0.0, 0.0]
      v0 = SA[v_perp, 0.0, 0.0]
      stateinit = [x0..., v0...]
      tspan = (0.0, 1.0)

      param = prepare(E_func, B_func; species = Ion(1, 1)) # Proton
      prob = ODEProblem(trace!, stateinit, tspan, param)
      sol = solve(prob, Tsit5())

      # Theoretical gyroradius
      # r = m * v_perp / (q * B)
      r_expected = m * v_perp / (q * B0)

      r_calc = get_gyroradius(sol, 0.5)
      @test r_calc ≈ r_expected rtol=1e-5
   end

   @testset "Relativistic" begin
      # Uniform B field
      B0 = 0.01
      B_func(x, t) = SA[0.0, 0.0, B0]
      E_func(x, t) = SA[0.0, 0.0, 0.0]

      # Particle parameters
      q = TestParticle.qₑ
      m = TestParticle.mₑ
      c = TestParticle.c

      # Initial state: relativistic velocity
      # Let v = 0.5c
      # γ = 1/sqrt(1-0.5^2) = 1/sqrt(0.75) ≈ 1.1547
      # u = γv
      v_mag = 0.5 * c
      gamma = 1.0 / sqrt(1 - (v_mag/c)^2)
      u_perp = gamma * v_mag

      x0 = SA[0.0, 0.0, 0.0]
      u0 = SA[u_perp, 0.0, 0.0] # γv
      stateinit = [x0..., u0...]
      tspan = (0.0, 1e-9)

      param = prepare(E_func, B_func; species = Electron)
      prob = ODEProblem(trace_relativistic!, stateinit, tspan, param)
      sol = solve(prob, Vern6())

      # Theoretical relativistic gyroradius
      # r = p_perp / (q * B) = m * (γv)_perp / (q * B)
      r_expected = m * u_perp / (abs(q) * B0)

      r_calc = get_gyroradius(sol, 0.5e-9)
      @test r_calc ≈ r_expected rtol=1e-5
   end
end
