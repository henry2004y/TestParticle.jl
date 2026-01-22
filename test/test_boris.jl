using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using OrdinaryDiffEq
using LinearAlgebra

@testset "Boris Solvers" begin
    # Definitions
    uniform_B2(x) = SA[0.0, 0.0, 0.01]
    function time_varying_B(x, t)
        Bz = if x[1] > 100t
            1.0e-8
        elseif x[1] < -10.0
            1.0e-8
        else
            0.0
        end
        return SA[0.0, 0.0, Bz]
    end
    zero_E = TP.ZeroField()

    function prob_func_boris_immutable(prob, i, repeat)
        # Note: prob.u0[5] = i*1e5 is not thread-safe!
        return prob = @views remake(prob; u0 = [prob.u0[1:4]..., i * 1.0e5, prob.u0[6]])
    end

    @testset "Basic Boris" begin
        x0 = [0.0, 0.0, 0.0]
        v0 = [0.0, 1.0e5, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 3.0e-8)
        dt = 3.0e-11
        param = prepare(zero_E, uniform_B2, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        sol = TP.solve(prob; dt, savestepinterval = 10)[1]

        @test sol.u[end] ≈ [
            -0.00010199139098074829, 3.4634030517007745e-5, 0.0,
            -60893.0154034644, -79322.38445151183, 0.0,
        ]
        @test length(sol.t) == length(sol.u)

        t = tspan[2] / 2
        @test sol(t) ≈ [
            -3.8587891411024776e-5, 5.3855910044312875e-5, 0.0,
            -94689.59405645168, 32154.016505320025, 0.0,
        ]

        prob = TraceProblem(stateinit, tspan, param; prob_func = prob_func_boris_immutable)
        trajectories = 4
        savestepinterval = 1000
        sols = TP.solve(prob, EnsembleThreads(); dt, savestepinterval, trajectories)
        @test sum(s -> sum(s.u[end][4]), sols) ≈ -608930.1540346438

        prob = TraceProblem(stateinit, tspan, param; prob_func = prob_func_boris_immutable)
        trajectories = 2
        savestepinterval = 1000
        sols = TP.solve(prob; dt, savestepinterval, trajectories)
        @test sum(s -> sum(s.u[end]), sols) ≈ -420646.1997670008

        x0 = [-1.0, 0.0, 0.0]
        v0 = [1.0e6, 0.0, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 0.01)
        dt = 1.0e-4
        param = prepare(zero_E, time_varying_B, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)
        sol = TP.solve(prob; dt, savestepinterval = 100)[1]
        @test sol[1, end] ≈ -512.8807058314515

        new_tspan = (0.0, 2.0e-8)
        new_prob = remake(prob; tspan = new_tspan)
        @test new_prob.tspan == new_tspan
    end

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

        @test sol_std[1].u[end][1] ≈ 10.0 atol = 1.0e-6
        @test sol_multi_2[1].u[end][1] ≈ 10.0 atol = 1.0e-6

        # Check velocities
        @test sol_std[1].u[end][4] ≈ 1.0 atol = 1.0e-6
        @test sol_multi_2[1].u[end][4] ≈ 1.0 atol = 1.0e-6

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
            @test hypot(v_end...) ≈ v0_norm atol = 2.0e-3
        end
    end

    @testset "Adaptive Boris" begin
        # Check constructor default
        alg1 = AdaptiveBoris(dtmax = 2.0)
        @test alg1.dtmin == 0.02

        # Check simulation
        x0 = [0.0, 0.0, 0.0]
        v0 = [0.0, 1.0e5, 0.0]
        stateinit = [x0..., v0...]
        # Gyroperiod for electron in 0.01T is ~3.57e-10 s
        # Let's run for 100 periods
        tperiod = abs(TP.get_gyroperiod(0.01; q = TP.qₑ, m = TP.mₑ))
        tspan = (0.0, 100 * tperiod)

        alg_adaptive = AdaptiveBoris(dtmax = tperiod, safety = 0.1)
        # dtmin should be 0.01 * tperiod

        param = prepare(zero_E, uniform_B2, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        sol = TP.solve(prob, alg_adaptive)[1]

        # Energy conservation check
        # |v| should be conserved
        v_end = sol.u[end][4:6]
        v_start = stateinit[4:6]
        @test norm(v_end) ≈ norm(v_start) atol = 1.0e-3 * norm(v_start)

        # Check that we took reasonable number of steps
        # Expect dt ~ safety / Omega = safety / (2pi/tperiod) = safety * tperiod / 2pi
        # dt ~ 0.1 * tperiod / 6.28 ~ 0.015 * tperiod
        # total / dt ~ 100 / 0.015 ~ 6666 steps
        # Warning: if length is 2, it means it triggered early exit or all steps rejected?
        @test 100 < length(sol.t) < 20000
    end
end
