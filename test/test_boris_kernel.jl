module test_boris_kernel

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using KernelAbstractions
using OrdinaryDiffEq
const KA = KernelAbstractions

@testset "Boris Kernel Solver" begin
    uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
    uniform_E(x) = SA[0.0, 0.0, 0.0]

    x0 = [0.0, 0.0, 0.0]
    v0 = [1.0e5, 0.0, 0.0]
    stateinit = [x0..., v0...]
    tspan = (0.0, 1.0e-6)
    dt = 1.0e-9

    param = prepare(uniform_E, uniform_B; species = Proton)
    prob = TraceProblem(stateinit, tspan, param)

    @testset "CPU Backend" begin
        backend = CPU()

        sol_gpu = TP.solve(prob, backend; dt, trajectories = 1, savestepinterval = 10)

        @test length(sol_gpu) == 1
        @test length(sol_gpu[1].t) > 0
        @test sol_gpu[1].t[1] == tspan[1]
        @test sol_gpu[1].t[end] ≈ tspan[2]

        x_final = sol_gpu[1].u[end][1]
        y_final = sol_gpu[1].u[end][2]

        @test abs(x_final) > 0
        @test abs(y_final) > 0
    end

    @testset "Multi-particle Kernel" begin
        backend = CPU()

        prob_func_gpu(prob, i, repeat) = remake(
            prob; u0 = [prob.u0[1:3]..., i * 1.0e4, 0.0, 0.0]
        )
        prob_multi = TraceProblem(stateinit, tspan, param; prob_func = prob_func_gpu)

        sols_gpu = TP.solve(prob_multi, backend; dt, trajectories = 5, savestepinterval = 100)

        @test length(sols_gpu) == 5

        for i in 1:5
            @test length(sols_gpu[i].t) > 0
            @test sols_gpu[i].t[1] == tspan[1]
        end
    end

    @testset "EnsembleThreads" begin
        backend = CPU()

        prob_func_gpu(prob, i, repeat) = remake(
            prob; u0 = [prob.u0[1:3]..., i * 1.0e4, 0.0, 0.0]
        )
        prob_multi = TraceProblem(stateinit, tspan, param; prob_func = prob_func_gpu)

        trajectories = 10
        sols_serial = TP.solve(
            prob_multi, backend, EnsembleSerial();
            dt, trajectories, savestepinterval = 100
        )
        sols_threads = TP.solve(
            prob_multi, backend, EnsembleThreads();
            dt, trajectories, savestepinterval = 100
        )

        @test length(sols_threads) == trajectories
        @test sum(sol.u[end] for sol in sols_threads) ≈ sum(sol.u[end] for sol in sols_serial)
    end

    @testset "Kernel vs Native Solver Equivalence" begin
        backend = CPU()

        sol_gpu = TP.solve(prob, backend; dt, trajectories = 1, savestepinterval = 10)
        sol_cpu = TP.solve(prob; dt, savestepinterval = 10)

        @test length(sol_gpu[1].t) == length(sol_cpu[1].t)

        # Check only first and last steps to avoid excessive test printing
        @test sol_gpu[1].t ≈ sol_cpu[1].t atol = 1.0e-6
        # Tight tolerances as the implementations should be algorithmically equivalent
        @test sol_gpu[1].u[end] ≈ sol_cpu[1].u[end] rtol = 1.0e-10 atol = 1.0e-10
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
        sol_gyro = TP.solve(prob_gyro, backend; dt = dt_gyro, savestepinterval = 10)

        # Check energy conservation at the final step to reduce test count
        vx = sol_gyro[1].u[end][4]
        vy = sol_gyro[1].u[end][5]
        vz = sol_gyro[1].u[end][6]
        v_mag = sqrt(vx^2 + vy^2 + vz^2)
        @test v_mag ≈ v_magnitude rtol = 1.0e-3
    end

    @testset "Saving Options" begin
        backend = CPU()

        sol_start_end = TP.solve(
            prob, backend; dt, trajectories = 1,
            save_everystep = false, save_start = true, save_end = true
        )
        @test length(sol_start_end[1].t) == 2
        @test sol_start_end[1].t[1] == tspan[1]
        @test sol_start_end[1].t[end] == tspan[2]

        sol_end_only = TP.solve(
            prob, backend; dt, trajectories = 1,
            save_everystep = false, save_start = false, save_end = true
        )
        @test length(sol_end_only[1].t) == 1
        @test sol_end_only[1].t[1] == tspan[2]
    end
end

end # module test_boris_kernel
