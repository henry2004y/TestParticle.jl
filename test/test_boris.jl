using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using OrdinaryDiffEq
using LinearAlgebra
using Distributed

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
    # Gradient B field: B = [0, 0, 0.01 * (1 + x)]
    # Strong gradient to test adaptive stepping
    gradient_B(x, t) = SA[0.0, 0.0, 0.01 * (1.0 + x[1])]

    # Constant E field for ExB drift
    # E = [0, 1e5, 0]
    constant_E(x, t) = SA[0.0, 1.0e5, 0.0]

    zero_E = ZeroField()

    constant_Ey(x, t) = SA[0.0, 1.0, 0.0]
    constant_Bz(x, t) = SA[0.0, 0.0, 1.0]

    function prob_func_boris_immutable(prob, ctx)
        # Note: prob.u0[5] = ctx.i*1e5 is not thread-safe!
        return prob = @views remake(prob; u0 = [prob.u0[1:4]..., ctx.i * 1.0e5, prob.u0[6]])
    end

    @testset "Basic Boris" begin
        x0 = [0.0, 0.0, 0.0]
        v0 = [0.0, 1.0e5, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 3.0e-8)
        dt = 3.0e-11
        param = prepare(zero_E, uniform_B2, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        sol = TP.solve(prob, Boris(); dt, savestepinterval = 10).u[1]

        @test sol.u[end] ≈ [
            -0.00010199137926394769, 3.46340469171306e-5, 0.0,
            -60893.043824907196, -79322.36263335795, 0.0,
        ]
        @test length(sol.t) == length(sol.u)

        t = tspan[2] / 2
        @test sol(t) ≈ [
            -3.8587882024745045e-5, 5.3855907133664956e-5, 0.0,
            -94689.58829601014, 32154.033469101283, 0.0,
        ]

        prob = TraceProblem(stateinit, tspan, param; prob_func = prob_func_boris_immutable)
        trajectories = 4
        savestepinterval = 1000
        sols = TP.solve(
            prob, Boris(), EnsembleThreads();
            dt, savestepinterval, trajectories
        )
        @test sum(s -> sum(s.u[end][4]), sols.u) ≈ -608930.4382490724

        prob = TraceProblem(stateinit, tspan, param; prob_func = prob_func_boris_immutable)
        trajectories = 2
        savestepinterval = 1000
        sols = TP.solve(
            prob, Boris(); dt, savestepinterval, trajectories
        )
        @test sum(s -> sum(s.u[end]), sols.u) ≈ -420646.2195768674

        x0 = [-1.0, 0.0, 0.0]
        v0 = [1.0e6, 0.0, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 0.01)
        dt = 1.0e-4
        param = prepare(zero_E, time_varying_B, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)
        sol = TP.solve(prob, Boris(); dt, savestepinterval = 100).u[1]
        @test sol[1, end] ≈ -512.8807214528281

        new_tspan = (0.0, 2.0e-8)
        new_prob = remake(prob; tspan = new_tspan)
        @test new_prob.tspan == new_tspan
    end

    @testset "Multistep Boris" begin
        # E cross B drift
        # E = (0, 1, 0), B = (0, 0, 1)
        # Analytic solution: particle moves at constant velocity v_drift = (1, 0, 0)
        # Position at t=10 should be (10, 0, 0)
        # q = 1, m = 1
        # The solver expects param to be (q2m, m, E, B, F)
        param = (1.0, 1.0, constant_Ey, constant_Bz, ZeroField())

        # Initial condition
        # Start at origin with drift velocity
        u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        tspan = (0.0, 10.0)
        dt = 0.1

        prob = TP.TraceProblem(u0, tspan, param)

        # Solve with standard Boris (n=1 by default)
        sol_std = TP.solve(prob, Boris(); dt)

        # 2-step
        sol_multi_2 = TP.solve(
            prob, MultistepBoris{2}(; n = 2); dt
        )

        @test sol_std.u[1].u[end][1] ≈ 10.0 atol = 1.0e-6
        @test sol_multi_2.u[1].u[end][1] ≈ 10.0 atol = 1.0e-6

        # Check velocities
        @test sol_std.u[1].u[end][4] ≈ 1.0 atol = 1.0e-6
        @test sol_multi_2.u[1].u[end][4] ≈ 1.0 atol = 1.0e-6

        # Test Gyrating particle
        # B = (0, 0, 1), E = 0
        # v = (1, 0, 0)
        # Gyroradius r = mv/qB = 1*1/1*1 = 1
        # Gyroperiod T = 2*pi*m/qB = 2*pi

        param_gyro = (1.0, 1.0, ZeroField(), constant_Bz, ZeroField())
        u0_gyro = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]

        prob_gyro = TP.TraceProblem(u0_gyro, (0.0, 2π), param_gyro)

        sol_1step_gyro = TP.solve(
            prob_gyro, Boris(); dt = 0.1
        )
        sol_2step_gyro = TP.solve(
            prob_gyro, MultistepBoris{2}(; n = 2); dt = 0.1
        )
        sol_4step_gyro = TP.solve(
            prob_gyro, MultistepBoris{2}(; n = 4); dt = 0.1
        )

        # After one period, should return to origin
        @test hypot(@views sol_1step_gyro.u[1].u[end][1:3]...) < 0.02
        @test hypot(@views sol_2step_gyro.u[1].u[end][1:3]...) < 0.02
        @test hypot(@views sol_4step_gyro.u[1].u[end][1:3]...) < 0.02

        # Energy conservation check (E=0 implies |v| is constant)
        v0_norm = hypot(@views u0_gyro[4:6]...)

        for sol in (sol_1step_gyro, sol_2step_gyro, sol_4step_gyro)
            v_end = @view sol.u[1].u[end][4:6]
            @test hypot(v_end...) ≈ v0_norm atol = 2.0e-3
        end
    end

    @testset "Hyper Boris" begin
        # E cross B drift tests
        param = (1.0, 1.0, constant_Ey, constant_Bz, ZeroField())
        u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        tspan = (0.0, 10.0)
        dt = 0.1

        prob = TP.TraceProblem(u0, tspan, param)

        # 2-step Hyper Boris (N=4, N=6)
        sol_hyper_4 = TP.solve(
            prob, MultistepBoris{4}(; n = 2); dt
        )
        sol_hyper_6 = TP.solve(
            prob, MultistepBoris{6}(; n = 2); dt
        )
        sol_hyper_single = TP.solve(
            prob, MultistepBoris{4}(; n = 1); dt
        )

        @test sol_hyper_4.u[1].u[end][1] ≈ 10.0 atol = 1.0e-6
        @test sol_hyper_6.u[1].u[end][1] ≈ 10.0 atol = 1.0e-6
        @test sol_hyper_single.u[1].u[end][1] ≈ 10.0 atol = 1.0e-6

        @test sol_hyper_4.u[1].u[end][4] ≈ 1.0 atol = 1.0e-6
        @test sol_hyper_6.u[1].u[end][4] ≈ 1.0 atol = 1.0e-6
        @test sol_hyper_single.u[1].u[end][4] ≈ 1.0 atol = 1.0e-6

        # Gyrating particle test
        param_gyro = (1.0, 1.0, ZeroField(), constant_Bz, ZeroField())
        u0_gyro = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        prob_gyro = TP.TraceProblem(u0_gyro, (0.0, 2π), param_gyro)

        sol_hyper_4_gyro = TP.solve(
            prob_gyro, MultistepBoris{4}(; n = 4); dt = 0.1
        )
        sol_hyper_6_gyro = TP.solve(
            prob_gyro, MultistepBoris{6}(; n = 4); dt = 0.1
        )

        @test hypot(@views sol_hyper_4_gyro.u[1].u[end][1:3]...) < 0.02
        @test hypot(@views sol_hyper_6_gyro.u[1].u[end][1:3]...) < 0.02

        v0_norm = hypot(@views u0_gyro[4:6]...)
        for sol in (sol_hyper_4_gyro, sol_hyper_6_gyro)
            v_end = @view sol.u[1].u[end][4:6]
            @test hypot(v_end...) ≈ v0_norm atol = 2.0e-3
        end
    end

    @testset "Adaptive Boris" begin
        # Check constructor default
        alg1 = AdaptiveBoris(safety = 0.1)
        @test alg1.safety == 0.1

        x0 = [0.0, 0.0, 0.0]
        v0 = [1.0e7, 0.0, 0.0]
        stateinit = [x0..., v0...]

        tperiod = abs(TP.get_gyroperiod(0.01; q = TP.qₑ, m = TP.mₑ))
        tspan = (0.0, 200 * tperiod)

        alg_adaptive = AdaptiveBoris(safety = 0.1)

        param = prepare(constant_E, gradient_B, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        sol = TP.solve(prob, alg_adaptive).u[1]

        dt_end = sol.t[end - 1] - sol.t[end - 2]
        dt_start = sol.t[11] - sol.t[10]

        @test dt_end < 0.8 * dt_start

        function total_energy(u)
            v = u[4:6]
            y = u[2]
            K = 0.5 * TP.mₑ * norm(v)^2 # Kinetic Energy
            U = TP.qₑ * (-1.0e5 * y) # Potential Energy
            return K + U
        end

        E_start = total_energy(sol.u[1])
        E_end = total_energy(sol.u[end])

        @test isapprox(E_end, E_start, rtol = 1.0e-3)
    end

    @testset "Nonzero tspan[1]" begin
        # B field that is only present when t > 5
        B_field(r, t) = t > 5 ? SA[0.0, 0.0, 0.01] : SA[0.0, 0.0, 0.0]
        E_field(r, t) = SA[0.0, 0.0, 0.0]

        param = prepare(E_field, B_field)

        # Start at t=10. If absolute time is used, B should be 0.01.
        tspan = (10.0, 10.1)
        dt = 0.001
        u0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
        prob = TraceProblem(u0, tspan, param)

        # Standard Boris
        sol_boris = TP.solve(prob, Boris(); dt = dt).u[1]
        # If B was 0.01, vx should have changed significantly
        @test abs(sol_boris.u[end][4]) < 1.0e5 - 100

        # Adaptive Boris
        alg_adaptive = AdaptiveBoris(safety = 0.1)
        sol_adaptive = TP.solve(prob, alg_adaptive).u[1]
        @test abs(sol_adaptive.u[end][4]) < 1.0e5 - 100
    end

    @testset "Output saving flags" begin
        # Setup
        x0 = [0.0, 0.0, 0.0]
        v0 = [0.0, 1.0e5, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 3.0e-8)
        dt = 3.0e-11
        zero_E = ZeroField()
        # uniform_B2(x) = SA[0.0, 0.0, 0.01] # Use global definition
        param = prepare(zero_E, uniform_B2, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        # Baseline: save_everystep=true (default), save_start=true (default implicit), save_end=true (default implicit)
        # savestepinterval = 10
        # nt = 1000. steps = 1000/10 = 100.
        # nout = 101 (0, 10, ..., 1000)
        sol = TP.solve(
            prob, Boris(); dt = dt, savestepinterval = 10
        ).u[1]
        @test length(sol.u) == 101
        @test sol.t[1] == tspan[1]
        @test sol.t[end] == tspan[2]

        # Scenario 2: Only final state
        sol = TP.solve(
            prob, Boris(); dt = dt, savestepinterval = 10,
            save_everystep = false, save_start = false,
            save_end = true
        ).u[1]
        @test length(sol.u) == 1
        @test sol.t[1] == tspan[2]

        # Scenario 3: Start and End
        sol = TP.solve(
            prob, Boris(); dt = dt, savestepinterval = 10,
            save_everystep = false, save_start = true,
            save_end = true
        ).u[1]
        @test length(sol.u) == 2
        @test sol.t[1] == tspan[1]
        @test sol.t[end] == tspan[2]

        # Scenario 4: Only start
        sol = TP.solve(
            prob, Boris(); dt = dt, savestepinterval = 10,
            save_everystep = false, save_start = true,
            save_end = false
        ).u[1]
        @test length(sol.u) == 1
        @test sol.t[1] == tspan[1]

        # Scenario 5: Every step but no start/end
        sol = TP.solve(
            prob, Boris(); dt = dt, savestepinterval = 10,
            save_everystep = true,
            save_start = false, save_end = false
        ).u[1]
        # Steps: 10, 20, ..., 990.
        @test length(sol.u) == 99
        @test sol.t[1] ≈ tspan[1] + 10 * dt
        @test sol.t[end] ≈ tspan[1] + 990 * dt

        # Scenario 6: Irregular interval
        sol = TP.solve(
            prob, Boris(); dt = dt, savestepinterval = 3,
            save_everystep = true, save_start = true,
            save_end = true
        ).u[1]
        @test length(sol.u) == 335
        @test sol.t[1] == tspan[1]
        @test sol.t[2] ≈ tspan[1] + 3 * dt
        @test sol.t[end - 1] ≈ tspan[1] + 999 * dt
        @test sol.t[end] == tspan[2]

        # Multistep Boris test
        sol_ms = TP.solve(
            prob, MultistepBoris{2}(; n = 2);
            dt = dt, savestepinterval = 10,
            save_everystep = false, save_start = true,
            save_end = true
        ).u[1]
        @test length(sol_ms.u) == 2
        @test sol_ms.t[1] == tspan[1]
        @test sol_ms.t[end] == tspan[2]

        # Save fields test
        sol_fields = TP.solve(
            prob, Boris();
            dt, savestepinterval = 10, save_fields = true
        ).u[1]
        # Check dimensions of the state vector (not the solution object length)
        @test length(sol_fields.u[1]) == 12
        # Check field values. E=0, B=[0,0,0.01].
        # Element 7-9 are E, 10-12 are B.
        @test sol_fields.u[1][7:9] == [0.0, 0.0, 0.0]
        @test sol_fields.u[1][10:12] == [0.0, 0.0, 0.01]
        @test sol_fields.u[end][7:9] == [0.0, 0.0, 0.0]
        @test sol_fields.u[end][10:12] == [0.0, 0.0, 0.01]
    end

    @testset "Save work" begin
        # Case 1: Gradient Drift Work
        # E = [0, 1e5, 0]
        # B = [0, 0, 0.01 * (1 + x)] => b=[0,0,1], ∇B=[0.01, 0, 0]
        # b x ∇B = [0, 0.01, 0]
        # P_grad ~ (b x ∇B) . E = 1000
        # P_par = 0 (E.b = 0)
        # P_fermi = 0 (b constant direction => κ=0)
        # P_betatron = 0 (static B)

        x0 = [0.0, 0.0, 0.0]
        v0 = [1.0e5, 0.0, 0.0] # v perp to B
        stateinit = [x0..., v0...]
        tspan = (0.0, 1.0e-9)
        dt = 1.0e-11

        param = prepare(constant_E, gradient_B, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        # Test save_work=true
        sol = TP.solve(
            prob, Boris();
            dt, savestepinterval = 1,
            save_work = true, save_everystep = true
        ).u[1]

        # Dimension check: 6 (state) + 4 (work) = 10
        @test length(sol.u[1]) == 10

        # Check values
        # Index 7: P_par
        # Index 8: P_fermi
        # Index 9: P_grad
        # Index 10: P_betatron

        work = sol.u[1][7:10]
        @test work[1] ≈ 0.0 atol = 1.0e-10 # P_par
        @test work[2] ≈ 0.0 atol = 1.0e-10 # P_fermi
        @test abs(work[3]) > 0.0       # P_grad should be non-zero
        @test work[4] ≈ 0.0 atol = 1.0e-10 # P_betatron

        # Case 2: Betatron Acceleration
        # B = [0, 0, 0.01 * (1 + 1e7 * t)]
        # E = 0
        # dB/dt = 0.01 * 1e7 = 1e5
        # P_betatron = mu * dB/dt

        function time_varying_B_linear(x, t)
            return SA[0.0, 0.0, 0.01 * (1.0 + 1.0e7 * t)]
        end

        param_beta = prepare(zero_E, time_varying_B_linear, species = Electron)
        prob_beta = TraceProblem(stateinit, tspan, param_beta)

        sol_beta = TP.solve(
            prob_beta, Boris();
            dt, savestepinterval = 1, save_work = true
        ).u[1]

        work_beta = sol_beta.u[1][7:10]
        @test work_beta[4] > 0.0 # P_betatron should be positive

        # Test with save_fields=true AND save_work=true
        sol_both = TP.solve(
            prob, Boris();
            dt, savestepinterval = 1,
            save_fields = true, save_work = true
        ).u[1]
        # Dim: 6 + 6 + 4 = 16
        @test length(sol_both.u[1]) == 16
        # E, B at 7-12
        # Work at 13-16
        @test sol_both.u[1][13:16] == sol.u[1][7:10]

        # Test Multistep Boris with save_work
        sol_ms = TP.solve(
            prob, MultistepBoris{2}(; n = 2);
            dt, save_work = true, savestepinterval = 1
        ).u[1]
        @test length(sol_ms.u[1]) == 10
        work_ms = sol_ms.u[1][7:10]
        @test work_ms[1] ≈ 0.0 atol = 1.0e-10
        @test abs(work_ms[3]) > 0.0

        # Test Adaptive Boris with save_work
        # Use simple AdaptiveBoris
        alg_adaptive = AdaptiveBoris(safety = 0.1)
        sol_adaptive = TP.solve(
            prob, alg_adaptive;
            save_work = true, save_everystep = true
        ).u[1]
        @test length(sol_adaptive.u[1]) == 10
        work_adaptive = sol_adaptive.u[1][7:10]
        @test work_adaptive[1] ≈ 0.0 atol = 1.0e-10
        @test abs(work_adaptive[3]) > 0.0
    end

    @testset "Post-processing fields and work" begin
        x0 = [0.0, 0.0, 0.0]
        v0 = [1.0e5, 0.0, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 1.0e-9)
        dt = 1.0e-11

        param = prepare(constant_E, gradient_B, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        # Reference solution with saved data
        sol_ref = TP.solve(
            prob, Boris();
            dt, savestepinterval = 1,
            save_fields = true, save_work = true
        ).u[1]

        # Solution without saved data
        sol = TP.solve(
            prob, Boris(); dt, savestepinterval = 1
        ).u[1]

        # Post-process
        E_post, B_post = get_fields(sol)
        work_post = get_work(sol)

        # Comparison
        # Reference indices:
        # 1-6: state
        # 7-9: E
        # 10-12: B
        # 13-16: work

        @test E_post[end] ≈ sol_ref.u[end][7:9]
        @test B_post[end] ≈ sol_ref.u[end][10:12]
        @test work_post[end] ≈ sol_ref.u[end][13:16]
    end

    @testset "Solver Limits" begin
        x0 = [0.0, 0.0, 0.0]
        v0 = [1.0e5, 0.0, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 1.0e-9)
        dt = 1.0e-11
        param = prepare(constant_E, gradient_B, species = Electron)
        prob = TraceProblem(stateinit, tspan, param)

        # maxiters limit (nt = 100)
        @test_throws ArgumentError TP.solve(
            prob, Boris(); dt, maxiters = 50
        )

        # min_dt limit
        dt_too_small = eps(Float64)
        @test_throws ArgumentError TP.solve(
            prob, Boris(); dt = dt_too_small
        )

        # Hyper Boris N parameter limit
        @test_throws ArgumentError MultistepBoris{1}()
        @test_throws ArgumentError MultistepBoris{3}(; n = 2)
    end

    @testset "EnsembleDistributed" begin
        pids = addprocs(2)
        try
            @everywhere pids using TestParticle
            @everywhere pids using StaticArrays

            x0 = [0.0, 0.0, 0.0]
            v0 = [0.0, 1.0e5, 0.0]
            stateinit = [x0..., v0...]
            tspan = (0.0, 3.0e-8)
            dt = 3.0e-11
            zero_E_dist = ZeroField()
            uniform_B2_dist(x) = SA[0.0, 0.0, 0.01]
            param_dist = prepare(zero_E_dist, uniform_B2_dist; species = Electron)

            function dist_prob_func(prob, ctx)
                return remake(
                    prob; u0 = SA[
                        prob.u0[1], prob.u0[2], prob.u0[3],
                        prob.u0[4], ctx.i * 1.0e5, prob.u0[6],
                    ]
                )
            end

            prob_dist = TraceProblem(
                stateinit, tspan, param_dist; prob_func = dist_prob_func
            )
            trajectories = 4
            savestepinterval = 1000

            sols_serial = TP.solve(
                prob_dist, Boris(), EnsembleSerial();
                dt, savestepinterval, trajectories
            )
            @testset "Boris" begin
                sols_dist = TP.solve(
                    prob_dist, Boris(),
                    EnsembleDistributed();
                    dt, savestepinterval, trajectories
                )

                @test length(sols_dist.u) == trajectories
                for i in 1:trajectories
                    @test sols_dist.u[i].u[end] ≈ sols_serial.u[i].u[end]
                end
            end
            @testset "SplitThreads" begin
                sols_split = TP.solve(
                    prob_dist, Boris(),
                    EnsembleSplitThreads();
                    dt, savestepinterval, trajectories
                )

                @test length(sols_split.u) == trajectories
                for i in 1:trajectories
                    @test sols_split.u[i].u[end] ≈ sols_serial.u[i].u[end]
                end
            end

            @testset "AdaptiveBoris" begin
                tperiod = abs(TP.get_gyroperiod(0.01; q = TP.qₑ, m = TP.mₑ))
                alg_adaptive = AdaptiveBoris(; safety = 0.1)

                sols_serial = TP.solve(
                    prob_dist, alg_adaptive, EnsembleSerial(); trajectories
                )

                @testset "EnsembleDistributed" begin
                    sols_dist = TP.solve(
                        prob_dist, alg_adaptive, EnsembleDistributed(); trajectories
                    )

                    @test length(sols_dist.u) == trajectories
                    for i in 1:trajectories
                        @test sols_dist.u[i].u[end] ≈ sols_serial.u[i].u[end]
                    end
                end

                @testset "EnsembleSplitThreads" begin
                    sols_split = TP.solve(
                        prob_dist, alg_adaptive, EnsembleSplitThreads(); trajectories
                    )

                    @test length(sols_split.u) == trajectories
                    for i in 1:trajectories
                        @test sols_split.u[i].u[end] ≈ sols_serial.u[i].u[end]
                    end
                end
            end
        finally
            rmprocs(pids)
        end
    end

end
