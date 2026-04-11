module test_isoutofdomain
using Test
using TestParticle
import TestParticle: solve, TraceProblem, TraceGCProblem, AdaptiveBoris, Proton, prepare
using StaticArrays
using LinearAlgebra
using OrdinaryDiffEq

@testset "isoutofdomain" begin

    @testset "NaN handling (Native Solvers)" begin
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

        # domain check: stop if r > 1.0
        isoutofdomain_r1(u, p, t) = norm(u[1:3]) > 1.0

        u0 = [0.9, 0.0, 0.0, 0.0, 1.0, 0.0] # Starts inside, moving out
        tspan = (0.0, 1.0)
        p = (1.0, 1.0, E_field_nan, B_field_nan, (x, t) -> 0.0)

        @testset "Boris" begin
            prob = TraceProblem(u0, tspan, p)
            sol = solve(prob, dt = 0.2, isoutofdomain = isoutofdomain_r1)

            # Check that no NaNs are in the solution
            for u in sol[1].u
                @test !any(isnan, u)
            end
            @test norm(sol[1].u[end][1:3]) <= 1.0
        end

        @testset "Adaptive Boris" begin
            prob = TraceProblem(u0, tspan, p)
            sol = solve(prob, AdaptiveBoris(safety = 0.1), isoutofdomain = isoutofdomain_r1)

            for u in sol[1].u
                @test !any(isnan, u)
            end
            @test norm(sol[1].u[end][1:3]) <= 1.0
        end

        @testset "GC RK4" begin
            u0_gc = [0.9, 0.0, 0.0, 0.0]
            p_gc = (1.0, 1.0, 0.0, E_field_nan, B_field_nan)
            prob = TraceGCProblem(u0_gc, tspan, p_gc)
            sol = solve(prob, dt = 0.2, alg = :rk4, isoutofdomain = isoutofdomain_r1)

            for u in sol[1].u
                @test !any(isnan, u)
            end
            @test norm(sol[1].u[end][1:3]) <= 1.0
        end

        @testset "GC RK45" begin
            u0_gc = [0.9, 0.0, 0.0, 0.0]
            p_gc = (1.0, 1.0, 0.0, E_field_nan, B_field_nan)
            prob = TraceGCProblem(u0_gc, tspan, p_gc)
            sol = solve(prob, dt = 0.01, alg = :rk45, isoutofdomain = isoutofdomain_r1)

            for u in sol[1].u
                @test !any(isnan, u)
            end
            @test norm(sol[1].u[end][1:3]) <= 1.0
        end
    end

    @testset "DiffEq Callback Integration" begin
        # Derived from runtests.jl
        isoutofdomain_r08(u, p, t) = hypot(u[1], u[2], u[3]) > 0.8

        x_grid = range(-10, 10, length = 15)
        y_grid = range(-10, 10, length = 20)
        z_grid = range(-10, 10, length = 25)
        B = fill(0.0, 3, length(x_grid), length(y_grid), length(z_grid))
        E = fill(0.0, 3, length(x_grid), length(y_grid), length(z_grid))

        B[3, :, :, :] .= 10.0e-9
        E[3, :, :, :] .= 5.0e-10

        # Match runtests.jl default (Proton)
        param = prepare(x_grid, y_grid, z_grid, E, B)

        u0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        tspan = (0.0, 1.0)

        prob = ODEProblem(trace!, u0, tspan, param)
        callback = DiscreteCallback(isoutofdomain_r08, terminate!)

        sol = OrdinaryDiffEq.solve(
            prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8,
            save_idxs = [1, 2, 3], callback, verbose = false
        )
        # Verify boundary termination
        # Tolerance is slightly higher for x86/ARM differences
        @test sol[1, end] ≈ 0.79411 rtol = 1.0e-2
        @test norm(sol.u[end][1:3]) <= 0.9
    end

    @testset "Spatial and Time Boundaries (Native Solvers)" begin
        # Derived from test_boris.jl and test_gc.jl
        B_field_fixed(r, t) = SA[1.0, 0.0, 0.0] # Parallel to v0 to ensure crossing
        E_field_fixed(r, t) = SA[0.0, 0.0, 0.0]
        param = prepare(E_field_fixed, B_field_fixed; species = Proton)

        @testset "Spatial rejection (AdaptiveBoris)" begin
            isoutofdomain_x(u, p, t) = u[1] > 0.5
            u0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0] # v = 1e5
            tspan = (0.0, 1.0e-5) # Should reach x=0.5 at t=5e-6
            prob = TraceProblem(u0, tspan, param)

            alg = AdaptiveBoris(safety = 0.1)
            sol = solve(prob, alg; isoutofdomain = isoutofdomain_x)[1]
            @test sol.u[end][1] <= 0.5 && sol.t[end] < tspan[2]
        end

        @testset "Time rejection (GC RK4)" begin
            # prepare_gc expects B(r) for time-independent fields
            B_func_gc(r) = SA[0.0, 0.0, 1.0]
            E_func_gc(r) = SA[0.0, 0.0, 0.0]

            u0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
            stateinit_gc, param_gc = TestParticle.prepare_gc(u0, E_func_gc, B_func_gc; species = Proton)
            tspan = (0.0, 1.0)
            prob = TraceGCProblem(stateinit_gc, tspan, param_gc)
            dt = 1.0e-4

            sol_early = solve(prob, dt = dt, alg = :rk4, isoutofdomain = (xv, p, t) -> t > 0.5)
            @test length(sol_early[1].t) == 5002
            @test sol_early[1].t[end] > 0.5
        end
    end

end

end # module test_isoutofdomain
