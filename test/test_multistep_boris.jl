import TestParticle as TP
using Test
using StaticArrays

@testset "Multistep Boris" begin
   # E cross B drift
   # E = (0, 1, 0), B = (0, 0, 1)
   # Analytic solution: particle moves at constant velocity v_drift = (1, 0, 0)
   # Position at t=10 should be (10, 0, 0)
   E(x, t) = SA[0.0, 1.0, 0.0]
   B(x, t) = SA[0.0, 0.0, 1.0]

   # q = 1, m = 1
   # The solver expects param to be (q2m, m, E, B)
   param = (1.0, 1.0, E, B)

   # Initial condition
   # Start at origin with drift velocity
   u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
   tspan = (0.0, 10.0)
   dt = 0.1

   prob = TP.TraceProblem(u0, tspan, param)

   # Solve with standard Boris (n=1 by default)
   sol_std = TP.solve(prob; dt)

   # 2-step
   sol_multi_2 = TP.solve(prob; dt, n = 2)

   @test sol_std[1].u[end][1]≈10.0 atol=1e-6
   @test sol_multi_2[1].u[end][1]≈10.0 atol=1e-6

   # Check velocities
   @test sol_std[1].u[end][4]≈1.0 atol=1e-6
   @test sol_multi_2[1].u[end][4]≈1.0 atol=1e-6

   # Test Gyrating particle
   # B = (0, 0, 1), E = 0
   # v = (1, 0, 0)
   # Gyroradius r = mv/qB = 1*1/1*1 = 1
   # Gyroperiod T = 2*pi*m/qB = 2*pi

   param_gyro = (1.0, 1.0, TP.ZeroField(), B)
   u0_gyro = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]

   prob_gyro = TP.TraceProblem(u0_gyro, (0.0, 2π), param_gyro)

   sol_1step_gyro = TP.solve(prob_gyro; dt = 0.1, n = 1)
   sol_2step_gyro = TP.solve(prob_gyro; dt = 0.1, n = 2)
   sol_4step_gyro = TP.solve(prob_gyro; dt = 0.1, n = 4)

   # After one period, should return to origin
   @test hypot(@views sol_1step_gyro[1].u[end][1:3]...) < 0.02
   @test hypot(@views sol_2step_gyro[1].u[end][1:3]...) < 0.02
   @test hypot(@views sol_4step_gyro[1].u[end][1:3]...) < 0.02

   # Energy conservation check (E=0 implies |v| is constant)
   v0_norm = hypot(@views u0_gyro[4:6]...)
   for sol in (sol_1step_gyro, sol_2step_gyro, sol_4step_gyro)
      v_end = @view sol[1].u[end][4:6]
      @test hypot(v_end...)≈v0_norm atol=1e-12
   end
end
