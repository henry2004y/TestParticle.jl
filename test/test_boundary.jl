module test_boundary

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using LinearAlgebra: norm
using OrdinaryDiffEq
using SciMLBase: ReturnCode

"Fields that turn to NaN outside the unit sphere."
E_nan(r, t) = norm(r) > 1.0 ? SVector{3}(NaN, NaN, NaN) : SVector{3}(0.0, 0.0, 0.0)
B_nan(r, t) = norm(r) > 1.0 ? SVector{3}(NaN, NaN, NaN) : SVector{3}(0.0, 0.0, 1.0)

"Stop once the state has gone NaN or left the unit sphere."
outside_or_nan(u, p, t) = any(isnan, u) || norm(u[1:3]) > 1.0

"Stop once the particle has passed x = 0.5."
past_half_x(u, p, t) = u[1] > 0.5

"Stop once the trace has passed t = 0.5."
past_half_time(u, p, t) = t > 0.5

"Uniform 1 T field along +x, for the spatial boundary checks."
B_along_x(r, t = 0.0) = SA[1.0, 0.0, 0.0]

"Uniform 1 T field along +z, for the guiding-center time check."
B_along_z(r) = SA[0.0, 0.0, 1.0]

"Zero electric field."
E_zero(r, t = 0.0) = SA[0.0, 0.0, 0.0]

@testset "Callbacks and Termination" begin

    @testset "NaN handling with DiscreteCallback" begin
        u0 = [0.9, 0.0, 0.0, 0.0, 1.0, 0.0] # Starts inside, moving out
        tspan = (0.0, 1.0)
        p = (1.0, 1.0, E_nan, B_nan, (x, t) -> 0.0)

        @testset "Boris" begin
            prob = TraceProblem(u0, tspan, p)
            sol = solve(
                prob, Boris(), dt = 0.2;
                isoutside = outside_or_nan
            )

            @test !any(isnan, sol.u[end])
            @test norm(sol.u[end][1:3]) <= 1.0
            @test sol.retcode == ReturnCode.Terminated
        end

        @testset "Adaptive Boris" begin
            prob = TraceProblem(u0, tspan, p)
            sol = solve(
                prob, AdaptiveBoris(safety = 0.05);
                isoutside = outside_or_nan
            )

            @test !any(isnan, sol.u[end])
            @test norm(sol.u[end][1:3]) <= 1.0
            @test sol.retcode == ReturnCode.Terminated
        end

        @testset "GC RK4" begin
            u0_gc = [0.9, 0.0, 0.0, 1.0]
            p_gc = (1.0, 1.0, 0.0, E_nan, B_nan)
            prob = TraceGCProblem(u0_gc, tspan, p_gc)
            sol = solve(prob, dt = 0.2, alg = :rk4; isoutside = outside_or_nan)

            @test !any(isnan, sol.u[1].u[end])
            @test norm(sol.u[1].u[end][1:3]) <= 1.0
            @test sol.u[1].retcode == ReturnCode.Terminated
        end

        @testset "GC RK45" begin
            u0_gc = [0.9, 0.0, 0.0, 1.0]
            p_gc = (1.0, 1.0, 0.0, E_nan, B_nan)
            prob = TraceGCProblem(u0_gc, tspan, p_gc)
            sol = solve(prob, dt = 0.01, alg = :rk45; isoutside = outside_or_nan)

            @test !any(isnan, sol.u[1].u[end])
            @test norm(sol.u[1].u[end][1:3]) <= 1.0
            @test sol.u[1].retcode == ReturnCode.Terminated
        end
    end

    @testset "Spatial and Time Boundaries (Native Solvers)" begin
        param = prepare(E_zero, B_along_x; species = Proton)

        @testset "Spatial rejection (AdaptiveBoris)" begin
            u0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
            tspan = (0.0, 1.0e-5)
            prob = TraceProblem(u0, tspan, param)

            alg = AdaptiveBoris(safety = 0.1)
            sol = solve(prob, alg; isoutside = past_half_x)
            @test sol.u[end][1] <= 0.5
            @test sol.t[end] < tspan[2]
            @test sol.retcode == ReturnCode.Terminated
        end

        @testset "Integer tspan (Boris)" begin
            u0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            tspan = (0, 10)
            prob = TraceProblem(u0, tspan, param)
            sol = solve(prob, Boris(); dt = 1.0)
            @test sol.t[end] == 10
            @test sol.retcode == ReturnCode.Success
        end

        @testset "Time rejection (GC RK4)" begin
            u0 = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
            stateinit_gc, param_gc = prepare_gc(u0, E_zero, B_along_z; species = Proton)
            tspan = (0.0, 1.0)
            prob = TraceGCProblem(stateinit_gc, tspan, param_gc)
            dt = 1.0e-4

            sol_early = solve(prob; dt, alg = :rk4, isoutside = past_half_time)
            @test sol_early.u[1].t[end] ≈ 0.5
            @test sol_early.u[1].retcode == ReturnCode.Terminated
        end
    end

    @testset "TerminateOutside convention" begin
        # The condition is written in TestParticle's `(u, p, t)` convention, so
        # the callback has to adapt it to SciML's `(u, t, integrator)`. Reading
        # `t` here pins that down: without the adaptation this `t` would be the
        # integrator instead, and the comparison would not be defined.
        callback = TerminateOutside(past_half_time)

        prob = ODEProblem((u, p, t) -> u, zeros(6), (0.0, 10.0), (1.0,))
        sol = solve(prob, Tsit5(); dt = 0.1, adaptive = false, callback)

        @test sol.retcode == ReturnCode.Terminated
        @test sol.t[end] ≈ 0.6
    end
end

end # module test_boundary
