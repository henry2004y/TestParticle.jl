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
end
