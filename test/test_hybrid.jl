using Test
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using SciMLBase

@testset "Hybrid Solver" begin

    # Define a simple analytic magnetic mirror field
    B0 = 1.0
    L = 2.0

    function B_simple(x, t)
        xx, yy, zz = x[1], x[2], x[3]
        Bx = -B0 * xx * zz / L^2
        By = -B0 * yy * zz / L^2
        Bz = B0 * (1 + zz^2 / L^2)
        return SVector(Bx, By, Bz)
    end

    # Also define 1-arg version for initialization in prepare_gc
    function B_simple(x)
        return B_simple(x, 0.0)
    end

    E_zero(x, t) = SVector(0.0, 0.0, 0.0)
    E_zero(x) = SVector(0.0, 0.0, 0.0)

    # Parameters
    m = TestParticle.Proton.m
    q = TestParticle.Proton.q
    Ek = 1.0e6 * 1.602e-19 # 1 MeV
    v0 = TestParticle.energy2velocity(Ek; m, q)

    # Initial condition
    x0 = SVector(0.5, 0.0, 0.5)

    # Velocity
    v_total = v0
    pitch_angle = 45.0
    v_perp = v_total * sind(pitch_angle)
    v_par = v_total * cosd(pitch_angle)

    v0_vec = SVector(0.0, v_perp, v_par)

    # Construct GC problem
    state_gc, params_gc = prepare_gc(vcat(x0, v0_vec), E_zero, B_simple; species = TestParticle.Proton)

    tspan = (0.0, 1.0e-4)
    prob_gc = ODEProblem(trace_gc!, state_gc, tspan, params_gc)

    @testset "Integration" begin
        # Test Hybrid
        sol = solve_hybrid(prob_gc, Tsit5(); epsilon = 0.1, dt = 1.0e-7)

        # Check return type
        @test sol isa SciMLBase.AbstractODESolution

        # Check content
        @test length(sol.t) > 0
        @test length(sol.u) == length(sol.t)
        @test sol.u[1] isa SVector{6, Float64}
    end

    @testset "Correctness" begin
        # 1. Compare with pure GC (large epsilon => always GC)
        sol_hybrid_gc = solve_hybrid(prob_gc, Tsit5(); epsilon = 1.0e10, dt = 1.0e-7)
        prob_pure_gc = ODEProblem(trace_gc!, state_gc, tspan, params_gc)
        sol_pure_gc = solve(prob_pure_gc, Tsit5(); dt = 1.0e-7)

        # Compare final positions (should be identical or very close)
        @test norm(sol_hybrid_gc.u[end][1:3] - sol_pure_gc.u[end][1:3]) < 1.0e-10

        # 2. Compare with pure Full (small epsilon => always Full)
        # Note: Hybrid starts in GC, checks condition, then switches.
        # If we set epsilon small, it should switch immediately.
        sol_hybrid_full = solve_hybrid(prob_gc, Tsit5(); epsilon = 1.0e-10, dt = 1.0e-7)

        # Construct equivalent pure full problem
        # We need to convert initial GC state to Full state
        state_full_init = TestParticle.gc_to_full(state_gc, params_gc, 0.0)
        # Construct full params: (q2m, m, E, B, F)
        p_full = (q / m, m, E_zero, B_simple, TestParticle.ZeroField())
        prob_pure_full = ODEProblem(trace, state_full_init, tspan, p_full)
        sol_pure_full = solve(prob_pure_full, Tsit5(); dt = 1.0e-7)

        # Compare final positions.
        # Relax tolerance slightly as one step might differ.
        @test norm(sol_hybrid_full.u[end][1:3] - sol_pure_full.u[end][1:3]) < 1.0e-4
    end

    @testset "Stability" begin
        # Test Energy Conservation during switching
        sol = solve_hybrid(prob_gc, Tsit5(); epsilon = 0.1, dt = 1.0e-7)


        # Check conservation of total energy if possible, or bounded behavior

        # Since we don't track mu output easily, we check stability
        @test all(x -> !isnan(x), Iterators.flatten(sol.u))
        @test Symbol(sol.retcode) in [:Success, :Terminated, :Default]
    end
end
