using TestParticle
using OrdinaryDiffEq
using StaticArrays
using Test
import TestParticle as TP

@testset "Field line tracing" begin
   @testset "Dipole field" begin
       param = prepare(TP.getE_dipole, TP.getB_dipole)
       L = 4.0 * TP.Rₑ
       stateinit = [L, 0.0, 0.0]

       # Check L shell conservation
       function check_L_shell(sol)
           errs = Float64[]
           for u in sol.u
               r, theta, phi = TP.cart2sph(u...)
               if abs(sin(theta)) > 1e-2 && r > TP.Rₑ
                  L_calc = r / sin(theta)^2
                  push!(errs, abs(L_calc - L))
               end
           end
           return isempty(errs) ? 0.0 : maximum(errs)
       end

       tspan = (0.0, 10.0 * TP.Rₑ)

       # Test helper function
       # Forward
       prob = trace_fieldline(stateinit, param, tspan; mode=:forward)
       sol = solve(prob, Tsit5())
       @test check_L_shell(sol) < 1e-3 * TP.Rₑ

       # Backward
       prob_back = trace_fieldline(stateinit, param, tspan; mode=:backward)
       sol_back = solve(prob_back, Tsit5())
       @test check_L_shell(sol_back) < 1e-3 * TP.Rₑ
       @test prob_back.tspan[2] < 0

       # Both
       probs = trace_fieldline(stateinit, param, tspan; mode=:both)
       @test length(probs) == 2
       sol1 = solve(probs[1], Tsit5())
       sol2 = solve(probs[2], Tsit5())
       @test check_L_shell(sol1) < 1e-3 * TP.Rₑ
       @test check_L_shell(sol2) < 1e-3 * TP.Rₑ

       # Test equation overloading (out-of-place equation)
       prob_op = ODEProblem(trace_fieldline, SA[stateinit...], tspan, param)
       sol_op = solve(prob_op, Tsit5())
       @test check_L_shell(sol_op) < 1e-3 * TP.Rₑ
   end

   @testset "Interpolated field" begin
      x = range(-10, 10, length=10)
      y = range(-10, 10, length=10)
      z = range(-10, 10, length=10)
      B = zeros(3, 10, 10, 10)
      B[1, :, :, :] .= 1.0 # Bx = 1

      param = prepare(x, y, z, zeros(3, 10, 10, 10), B)
      stateinit = [0.0, 0.0, 0.0]
      tspan = (0.0, 5.0)

      prob = trace_fieldline(stateinit, param, tspan; mode=:forward)
      sol = solve(prob, Tsit5())

      @test sol.u[end][1] ≈ 5.0 atol=1e-2
      @test sol.u[end][2] ≈ 0.0 atol=1e-6
      @test sol.u[end][3] ≈ 0.0 atol=1e-6
   end
end
