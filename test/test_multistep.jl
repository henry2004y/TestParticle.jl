using TestParticle
using Test
using StaticArrays
using LinearAlgebra

@testset "Multistep Boris" begin
    # E cross B drift
    # E = (0, 1, 0), B = (0, 0, 1) -> v_drift = (1, 0, 0)

    E(x, t) = SA[0.0, 1.0, 0.0]
    B(x, t) = SA[0.0, 0.0, 1.0]

    # q = 1, m = 1
    # The solver expects param to be (q2m, m, E, B) or similar structure where it extracts 1, 3, 4.
    # Actually update_velocity! does: q2m, _, Efunc, Bfunc = param
    param = (1.0, 1.0, E, B)

    # Initial condition
    # Start at origin with drift velocity
    u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    tspan = (0.0, 10.0)
    dt = 0.1

    prob = TraceProblem(u0, tspan, param)

    # Solve with standard Boris (n=1 by default)
    sol_std = TestParticle.solve(prob; dt=dt, n=1)

    # Solve with Multistep Boris (n=1) explicitly
    # Since I modified solve to take n, and n=1 dispatches to standard boris,
    # I need to force use of multistep solver to test it for n=1.
    # But I made n=1 call _boris!. So to test _multistep_boris! with n=1,
    # I need to modify pusher.jl temporarily or access the internal function.
    # However, for n > 1, it uses multistep.

    # Let's test n=2
    sol_multi_2 = TestParticle.solve(prob; dt=dt, n=2)

    # Analytic solution: particle moves at constant velocity (1, 0, 0)
    # Position at t=10 should be (10, 0, 0)

    @test sol_std[1].u[end][1] ≈ 10.0 atol=1e-2
    @test sol_multi_2[1].u[end][1] ≈ 10.0 atol=1e-2

    # Check velocities
    @test sol_std[1].u[end][4] ≈ 1.0 atol=1e-4
    @test sol_multi_2[1].u[end][4] ≈ 1.0 atol=1e-4

    # Test Gyrating particle
    # B = (0, 0, 1), E = 0
    # v = (1, 0, 0)
    # Gyroradius r = mv/qB = 1*1/1*1 = 1
    # Gyroperiod T = 2*pi*m/qB = 2*pi

    E_zero(x, t) = SA[0.0, 0.0, 0.0]
    param_gyro = (1.0, 1.0, E_zero, B)
    u0_gyro = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0] # start at origin, v=(1,0,0)
    # This corresponds to gyration center at (0, 1, 0)?
    # Force F = q v x B = 1 * (1,0,0) x (0,0,1) = (0, -1, 0).
    # Acceleration a = (0, -1, 0).
    # Trajectory should be a circle tangent to x-axis at origin, center at (0, 1, 0)? No.
    # v x B = (vy*Bz - vz*By, vz*Bx - vx*Bz, vx*By - vy*Bx)
    #       = (0 - 0, 0 - 1*1, 0) = (0, -1, 0)
    # Acceleration is downwards.
    # If center is (0, 1, 0), radius 1.
    # x = r sin(wt), y = r (1 - cos(wt)) ?
    # at t=0, x=0, y=0. v=(r w, 0) = (1, 0). Correct.

    prob_gyro = TraceProblem(u0_gyro, (0.0, 2π), param_gyro)

    sol_std_gyro = TestParticle.solve(prob_gyro; dt=0.1, n=1)
    sol_multi_2_gyro = TestParticle.solve(prob_gyro; dt=0.1, n=2)
    sol_multi_4_gyro = TestParticle.solve(prob_gyro; dt=0.1, n=4)

    # After one period, should return to origin (approx)
    @test norm(sol_std_gyro[1].u[end][1:3]) < 0.1
    @test norm(sol_multi_2_gyro[1].u[end][1:3]) < 0.1
    @test norm(sol_multi_4_gyro[1].u[end][1:3]) < 0.1

    # Energy conservation check (E=0 implies |v| is constant)
    v0_norm = norm(u0_gyro[4:6])
    for sol in [sol_std_gyro, sol_multi_2_gyro, sol_multi_4_gyro]
        v_end = sol[1].u[end][4:6]
        @test norm(v_end) ≈ v0_norm atol=1e-12
    end

    # Verify that Multistep n=1 is consistent with Standard Boris
    # Note: solve with n=1 dispatches to _boris!, so we can't easily test _multistep_boris! with n=1 via solve.
    # But we can call update_velocity_multistep! manually to check against single step update.

    paramBoris = TestParticle.MultistepBorisMethod()
    xv = MVector{6, Float64}(u0_gyro)
    dt_step = 0.1
    t = 0.0
    # Use n=1
    TestParticle.update_velocity_multistep!(xv, paramBoris, param_gyro, dt_step, t, 1)
    v_multi_1 = copy(xv[4:6])

    # Standard Boris update
    xv_std = MVector{6, Float64}(u0_gyro)
    paramBorisStd = TestParticle.BorisMethod()
    TestParticle.update_velocity!(xv_std, paramBorisStd, param_gyro, dt_step, t)
    v_std = copy(xv_std[4:6])

    @test v_multi_1 ≈ v_std atol=1e-14

end
