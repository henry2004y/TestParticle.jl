using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using LinearAlgebra: norm
using OrdinaryDiffEq: ReturnCode

@testset "Callbacks and Termination" begin

    @testset "NaN handling with DiscreteCallback" begin
        # Define a field that returns NaN outside R < 1.0
        function E_field_nan(r, t)
            if norm(r) > 1.0
                return SVector{3}(NaN, NaN, NaN)
            else
                return SVector{3}(0.0, 0.0, 0.0)
            end
        end

        function B_field_nan(r, t)
            if norm(r) > 1.0
                return SVector{3}(NaN, NaN, NaN)
            else
                return SVector{3}(0.0, 0.0, 1.0)
            end
        end

        # domain check: stop if r > 1.0 or if NaN is encountered
        isoutside(u, p, t) = any(isnan, u) || norm(u[1:3]) > 1.0
        callback = TerminateOutside(isoutside)

        u0 = [0.9, 0.0, 0.0, 0.0, 1.0, 0.0] # Starts inside, moving out
        tspan = (0.0, 1.0)
        p = (1.0, 1.0, E_field_nan, B_field_nan, (x, t) -> 0.0)

        @testset "Boris" begin
            prob = TraceProblem(u0, tspan, p)
            sol = TP.solve(prob, dt = 0.2; isoutside = callback.condition)[1]

            @test !any(isnan, sol.u[end])
            @test norm(sol.u[end][1:3]) <= 1.0
            @test sol.retcode == ReturnCode.Terminated
        end

        @testset "Adaptive Boris" begin
            prob = TraceProblem(u0, tspan, p)
            sol = TP.solve(
                prob, AdaptiveBoris(safety = 0.05);
                isoutside = callback.condition
            )[1]

            @test !any(isnan, sol.u[end])
            @test norm(sol.u[end][1:3]) <= 1.0
            @test sol.retcode == ReturnCode.Terminated
        end

        @testset "GC RK4" begin
            u0_gc = [0.9, 0.0, 0.0, 1.0]
            p_gc = (1.0, 1.0, 0.0, E_field_nan, B_field_nan)
            prob = TraceGCProblem(u0_gc, tspan, p_gc)
            sol = TP.solve(prob, dt = 0.2, alg = :rk4; isoutside = callback.condition)[1]

            @test !any(isnan, sol.u[end])
            @test norm(sol.u[end][1:3]) <= 1.0
            @test sol.retcode == ReturnCode.Terminated
        end

        @testset "GC RK45" begin
            u0_gc = [0.9, 0.0, 0.0, 1.0]
            p_gc = (1.0, 1.0, 0.0, E_field_nan, B_field_nan)
            prob = TraceGCProblem(u0_gc, tspan, p_gc)
            sol = TP.solve(prob, dt = 0.01, alg = :rk45; isoutside = callback.condition)[1]

            @test !any(isnan, sol.u[end])
            @test norm(sol.u[end][1:3]) <= 1.0
            @test sol.retcode == ReturnCode.Terminated
        end
    end

    @testset "Spatial and Time Boundaries (Native Solvers)" begin
        B_field_fixed(r, t = 0.0) = SA[1.0, 0.0, 0.0]
        E_field_fixed(r, t = 0.0) = SA[0.0, 0.0, 0.0]
        param = prepare(E_field_fixed, B_field_fixed; species = Proton)

        @testset "Spatial rejection (AdaptiveBoris)" begin
            isoutside(u, p, t) = u[1] > 0.5
            u0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
            tspan = (0.0, 1.0e-5)
            prob = TraceProblem(u0, tspan, param)

            alg = AdaptiveBoris(safety = 0.1)
            sol = TP.solve(prob, alg; isoutside)[1]
            @test sol.u[end][1] <= 0.5
            @test sol.t[end] < tspan[2]
            @test sol.retcode == ReturnCode.Terminated
        end

        @testset "Integer tspan (Boris)" begin
            u0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            tspan = (0, 10)
            prob = TraceProblem(u0, tspan, param)
            sol = TP.solve(prob; dt = 1.0)[1]
            @test sol.t[end] == 10
            @test sol.retcode == ReturnCode.Success
        end

        @testset "Time rejection (GC RK4)" begin
            B_func_gc(r) = SA[0.0, 0.0, 1.0]
            E_func_gc(r) = SA[0.0, 0.0, 0.0]

            u0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
            stateinit_gc, param_gc = prepare_gc(u0, E_func_gc, B_func_gc; species = Proton)
            tspan = (0.0, 1.0)
            prob = TraceGCProblem(stateinit_gc, tspan, param_gc)
            dt = 1.0e-4

            isoutside = (u, p, t) -> t > 0.5
            sol_early = TP.solve(prob; dt, alg = :rk4, isoutside)[1]
            @test sol_early.t[end] ≈ 0.5
            @test sol_early.retcode == ReturnCode.Terminated
        end
    end
end
