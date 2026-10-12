module test_boris_kernel

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using KernelAbstractions
using OrdinaryDiffEq
using SciMLBase: ReturnCode
using ..test_common: uniform_B, zero_E
const KA = KernelAbstractions

"A distinct initial speed per particle, so reordering would be visible."
prob_func_gpu(prob, ctx) = remake(
    prob; u0 = [prob.u0[1:3]..., ctx.sim_id * 1.0e4, 0.0, 0.0]
)

@testset "Boris Kernel Solver" begin
    x0 = [0.0, 0.0, 0.0]
    v0 = [1.0e5, 0.0, 0.0]
    stateinit = [x0..., v0...]
    tspan = (0.0, 1.0e-6)
    dt = 1.0e-9

    param = prepare(zero_E, uniform_B; species = Proton)
    prob = TraceProblem(stateinit, tspan, param)

    @testset "CPU Backend" begin
        backend = CPU()

        sol_gpu = solve(prob, Boris(), backend; dt, trajectories = 1, saveat = 10 * dt)

        @test length(sol_gpu.u) == 1
        @test length(sol_gpu.u[1].t) > 0
        @test sol_gpu.u[1].t[1] == tspan[1]
        @test sol_gpu.u[1].t[end] ≈ tspan[2]

        x_final = sol_gpu.u[1].u[end][1]
        y_final = sol_gpu.u[1].u[end][2]

        @test abs(x_final) > 0
        @test abs(y_final) > 0
    end

    @testset "Multi-particle Kernel" begin
        backend = CPU()

        prob_multi = TraceProblem(stateinit, tspan, param; prob_func = prob_func_gpu)

        sols_gpu = solve(prob_multi, Boris(), backend; dt, trajectories = 5, saveat = 100 * dt)

        @test length(sols_gpu.u) == 5

        for i in 1:5
            @test length(sols_gpu.u[i].t) > 0
            @test sols_gpu.u[i].t[1] == tspan[1]
        end
    end

    @testset "EnsembleThreads" begin
        backend = CPU()

        prob_multi = TraceProblem(stateinit, tspan, param; prob_func = prob_func_gpu)

        trajectories = 10
        sols_serial = solve(
            prob_multi, Boris(), backend, EnsembleSerial();
            dt, trajectories, saveat = 100 * dt
        )
        sols_threads = solve(
            prob_multi, Boris(), backend, EnsembleThreads();
            dt, trajectories, saveat = 100 * dt
        )

        @test length(sols_threads.u) == trajectories
        @test sum(sol.u[end] for sol in sols_threads.u) ≈ sum(sol.u[end] for sol in sols_serial.u)
    end

    @testset "Kernel vs Native Solver Equivalence" begin
        backend = CPU()

        sol_gpu = solve(prob, Boris(), backend; dt, trajectories = 1, saveat = 10 * dt)
        sol_cpu = solve(prob, Boris(); dt, saveat = 10 * dt)

        @test length(sol_gpu.u[1].t) == length(sol_cpu.t)

        # Check only first and last steps to avoid excessive test printing
        @test sol_gpu.u[1].t ≈ sol_cpu.t atol = 1.0e-6
        # Tight tolerances as the implementations should be algorithmically equivalent
        @test sol_gpu.u[1].u[end] ≈ sol_cpu.u[end] rtol = 1.0e-10 atol = 1.0e-10
    end

    @testset "Energy Conservation" begin
        backend = CPU()

        B_magnitude = 1.0e-8
        v_magnitude = 1.0e5

        gyrofrequency = abs(TP.qᵢ) * B_magnitude / TP.mᵢ
        gyroperiod = 2π / gyrofrequency

        tspan_gyro = (0.0, 2 * gyroperiod)
        dt_gyro = gyroperiod / 100

        prob_gyro = TraceProblem(stateinit, tspan_gyro, param)
        sol_gyro = solve(prob_gyro, Boris(), backend; dt = dt_gyro, saveat = 10 * dt_gyro)

        # Check energy conservation at the final step to reduce test count
        vx = sol_gyro.u[1].u[end][4]
        vy = sol_gyro.u[1].u[end][5]
        vz = sol_gyro.u[1].u[end][6]
        v_mag = sqrt(vx^2 + vy^2 + vz^2)
        @test v_mag ≈ v_magnitude rtol = 1.0e-3
    end

    @testset "Saving Options" begin
        backend = CPU()

        sol_start_end = solve(
            prob, Boris(), backend; dt, trajectories = 1,
            save_everystep = false, save_start = true, save_end = true
        )
        @test length(sol_start_end.u[1].t) == 2
        @test sol_start_end.u[1].t[1] == tspan[1]
        @test sol_start_end.u[1].t[end] == tspan[2]

        sol_end_only = solve(
            prob, Boris(), backend; dt, trajectories = 1,
            save_everystep = false, save_start = false, save_end = true
        )
        @test length(sol_end_only.u[1].t) == 1
        @test sol_end_only.u[1].t[1] == tspan[2]
    end

    @testset "Solver Limits" begin
        backend = CPU()

        # maxiters limit (nt = 1000)
        @test_throws ArgumentError solve(prob, Boris(), backend; dt, maxiters = 500)

        # min_dt limit
        dt_too_small = eps(Float64)
        @test_throws ArgumentError solve(prob, Boris(), backend; dt = dt_too_small)
    end

    @testset "Saveat" begin
        backend = CPU()

        sol_all = solve(prob, Boris(), backend; dt, trajectories = 1).u[1]

        # Requesting intermediate output must not change the integration, only
        # where the state is reported. The requested times land on step
        # boundaries here, so the interpolant has to reproduce them exactly.
        ts = collect(1.0e-7:1.0e-7:9.0e-7)
        sol_at = solve(
            prob, Boris(), backend; dt, trajectories = 1, saveat = ts
        ).u[1]

        @test sol_at.t ≈ vcat(0.0, ts, 1.0e-6)
        for (k, t) in enumerate(sol_at.t)
            @test sol_at.u[k] ≈ sol_all(t)
        end

        # An interval is accepted too: only its interior is used, since the ends
        # of the span are reported separately.
        sol_interval = solve(
            prob, Boris(), backend; dt, trajectories = 1, saveat = 2.5e-7
        ).u[1]
        @test sol_interval.t ≈ collect(0.0:2.5e-7:1.0e-6)

    end

    @testset "Out-of-domain Termination" begin
        backend = CPU()

        # With both fields zero the state moves in a straight line, so x = v t.
        # The boundary x = 0.505 falls between the steps ending at 5.0e-6 and
        # 5.1e-6, so the step ending at 5.1e-6 is the one that leaves the domain.
        stateinit = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
        tspan = (0.0, 1.0e-5)
        dt = 1.0e-7

        param = prepare(zero_E, zero_E; species = Proton)
        prob = TraceProblem(stateinit, tspan, param)

        "Stop once the particle has passed x = 0.505."
        past_half_x(u, p, t) = u[1] > 0.505

        @testset "every step" begin
            sol = solve(prob, Boris(), backend; dt, isoutside = past_half_x).u[1]

            @test sol.retcode == ReturnCode.Terminated
            # The step that left the domain is rolled back, so the trace ends at
            # the last state inside it.
            @test sol.t[end] ≈ 5.0e-6
            @test sol.t[end] < tspan[2]
            @test sol.u[end][1] ≈ 0.5
            @test all(u -> u[1] <= 0.505, sol.u)
        end

        @testset "saveat" begin
            sol = solve(
                prob, Boris(), backend; dt, isoutside = past_half_x, saveat = 3.0e-7
            ).u[1]

            @test sol.retcode == ReturnCode.Terminated
            @test sol.t ≈ vcat(0.0, 3.0e-7:3.0e-7:4.8e-6, 5.0e-6)
            @test sol.u[end][1] ≈ 0.5
            @test all(u -> u[1] <= 0.505, sol.u)
        end

        @testset "end states only" begin
            sol = solve(
                prob, Boris(), backend; dt, isoutside = past_half_x,
                save_everystep = false
            ).u[1]

            @test sol.retcode == ReturnCode.Terminated
            @test sol.t ≈ [0.0, 5.0e-6]
            @test sol.u[end][1] ≈ 0.5
        end

        @testset "a trace that stays inside runs to completion" begin
            sol = solve(
                prob, Boris(), backend; dt, isoutside = (u, p, t) -> false
            ).u[1]

            @test sol.retcode == ReturnCode.Success
            @test sol.t[end] ≈ tspan[2]
            @test sol.u[end][1] ≈ 1.0
        end

        @testset "the check does not change a trace that stays inside" begin
            sol_plain = solve(prob, Boris(), backend; dt).u[1]
            sol_checked = solve(
                prob, Boris(), backend; dt, isoutside = (u, p, t) -> false
            ).u[1]

            @test sol_checked.t ≈ sol_plain.t
            @test sol_checked.u[end] ≈ sol_plain.u[end]
        end

        @testset "raw output" begin
            raw = solve(
                prob, Boris(), backend; dt, isoutside = past_half_x, raw_output = true
            )

            @test raw.terminated == [true]
            # 51 states are saved, the start plus the 50 steps that stayed inside.
            @test raw.u[1, 1, 51] ≈ 0.5
            # The slots the trace never reached hold its last state.
            @test raw.u[1, 1, 52] ≈ 0.5
            @test raw.u[1, 1, end] ≈ 0.5

            raw_ok = solve(prob, Boris(), backend; dt, raw_output = true)
            @test raw_ok.terminated == [false]
            @test raw_ok.u[1, 1, end] ≈ 1.0
        end

        @testset "ensemble with mixed outcomes" begin
            # The speed grows with the trajectory index, so only the faster ones
            # reach the boundary within the time span.
            prob_func = (p, ctx) -> remake(
                p; u0 = [0.0, 0.0, 0.0, ctx.sim_id * 1.0e4, 0.0, 0.0]
            )
            prob_multi = TraceProblem(stateinit, tspan, param; prob_func)

            sols = solve(
                prob_multi, Boris(), backend; dt, isoutside = past_half_x, trajectories = 8
            )

            @test [s.retcode for s in sols.u] == [
                ReturnCode.Success, ReturnCode.Success, ReturnCode.Success,
                ReturnCode.Success, ReturnCode.Success, ReturnCode.Terminated,
                ReturnCode.Terminated, ReturnCode.Terminated,
            ]
            for s in sols.u
                if s.retcode == ReturnCode.Terminated
                    @test s.t[end] < tspan[2]
                    @test s.u[end][1] <= 0.505
                end
            end

            raw = solve(
                prob_multi, Boris(), backend; dt, isoutside = past_half_x,
                trajectories = 8, raw_output = true
            )
            @test raw.terminated ==
                [false, false, false, false, false, true, true, true]
            for i in 6:8
                @test raw.u[i, 1, end] <= 0.505
                @test raw.u[i, 1, end] == raw.u[i, 1, end - 1]
            end
        end

        @testset "a trace that starts outside stops at once" begin
            prob_out = TraceProblem([2.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan, param)
            sol = solve(
                prob_out, Boris(), backend; dt, isoutside = past_half_x,
                save_everystep = false
            ).u[1]

            @test sol.retcode == ReturnCode.Terminated
            @test sol.u[end][1] == 2.0
            # The start state is the only one inside the domain, so the end slot
            # would repeat it and is dropped instead.
            @test length(sol.t) == 1
            @test sol.t == [tspan[1]]
        end

        @testset "agreement with the native solver" begin
            sol_kernel = solve(prob, Boris(), backend; dt, isoutside = past_half_x).u[1]
            sol_native = solve(prob, Boris(); dt, isoutside = past_half_x)

            # Both stop at the last state inside the domain. The native loop also
            # saves the step that left it before rolling that step back, so it
            # reports two more states, the last of them the same one.
            @test sol_kernel.u[end] ≈ sol_native.u[end]
            @test length(sol_native.t) == length(sol_kernel.t) + 2
            @test sol_native.retcode == ReturnCode.Terminated
        end

        @testset "backward time span" begin
            # The state moves toward +x as t runs backward, so the boundary is
            # crossed at t = -5.0e-6.
            prob_bw = TraceProblem(
                [0.0, 0.0, 0.0, -1.0e5, 0.0, 0.0], (0.0, -1.0e-5), param
            )
            sol = solve(
                prob_bw, Boris(), backend; dt = -1.0e-7, isoutside = past_half_x
            ).u[1]

            @test sol.retcode == ReturnCode.Terminated
            @test sol.t[end] ≈ -5.0e-6
            @test sol.t[end] > -1.0e-5
            @test sol.u[end][1] ≈ 0.5
            @test all(u -> u[1] <= 0.505, sol.u)
        end
    end
end

end # module test_boris_kernel
