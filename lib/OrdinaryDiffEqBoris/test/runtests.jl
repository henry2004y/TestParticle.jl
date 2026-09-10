using Test
using OrdinaryDiffEqBoris
using StaticArrays
using LinearAlgebra: norm

struct CustomParam{E, B}
    q2m::Float64
    m::Float64
    E::E
    B::B
end

OrdinaryDiffEqBoris.get_q2m(p::CustomParam) = p.q2m
OrdinaryDiffEqBoris.get_EField(p::CustomParam) = p.E
OrdinaryDiffEqBoris.get_BField(p::CustomParam) = p.B

@testset "OrdinaryDiffEqBoris.jl" begin
    # Definitions
    constant_Ey(x, t) = SA[0.0, 1.0, 0.0]
    constant_Bz(x, t) = SA[0.0, 0.0, 1.0]
    ZeroField() = (x, t) -> SA[0.0, 0.0, 0.0]

    dummy_f(u, p, t) = u # ODEProblem requires an f, but Boris uses p directly

    @testset "E cross B drift - Standard Boris" begin
        # E = (0, 1, 0), B = (0, 0, 1)
        # Analytic solution: particle moves at constant velocity v_drift = (1, 0, 0)
        param = (1.0, 1.0, constant_Ey, constant_Bz, ZeroField())
        u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        tspan = (0.0, 10.0)
        prob = ODEProblem(dummy_f, u0, tspan, param)

        dt = 0.1
        # Test standard Boris
        sol_boris = solve(prob, Boris(); dt)
        @test sol_boris.u[end][1] ≈ 10.0 atol = 1.0e-6
        @test sol_boris.u[end][4] ≈ 1.0 atol = 1.0e-6
    end

    @testset "E cross B drift - Multistep Boris" begin
        param = (1.0, 1.0, constant_Ey, constant_Bz, ZeroField())
        u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        tspan = (0.0, 10.0)
        prob = ODEProblem(dummy_f, u0, tspan, param)

        dt = 0.1
        sol_multi_2 = solve(prob, MultistepBoris2(n = 2); dt)
        @test sol_multi_2.u[end][1] ≈ 10.0 atol = 1.0e-6
        @test sol_multi_2.u[end][4] ≈ 1.0 atol = 1.0e-6

        # Test Hyper Boris N=4
        sol_hyper_4 = solve(prob, MultistepBoris4(n = 2); dt)
        @test sol_hyper_4.u[end][1] ≈ 10.0 atol = 1.0e-6
    end

    @testset "Gyrating particle" begin
        param_gyro = (1.0, 1.0, ZeroField(), constant_Bz, ZeroField())
        u0_gyro = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        prob_gyro = ODEProblem(dummy_f, u0_gyro, (0.0, 2π), param_gyro)

        # After one period, should return to origin
        sol_1step_gyro = solve(prob_gyro, MultistepBoris2(n = 1); dt = 0.1)
        @test hypot(@views sol_1step_gyro.u[end][1:3]...) < 0.02

        sol_hyper_6_gyro = solve(prob_gyro, MultistepBoris6(n = 4); dt = 0.1)
        @test hypot(@views sol_hyper_6_gyro.u[end][1:3]...) < 0.02
    end

    @testset "Adaptive Boris" begin
        # The electron-like setup below keeps the E cross B drift comparable to
        # the particle speed. In a drift-dominated regime the kinetic and
        # potential energies nearly cancel, which amplifies the energy error and
        # makes a relative tolerance on the total energy meaningless.
        q2m = -1.75882001076e11
        m = 9.1093837015e-31
        constant_E(x, t) = SA[0.0, 1.0e5, 0.0]
        gradient_B(x, t) = SA[0.0, 0.0, 0.01 * (1.0 + x[1])]
        param = (q2m, m, constant_E, gradient_B, ZeroField())

        u0 = SA[0.0, 0.0, 0.0, 1.0e7, 0.0, 0.0]
        safety = 0.1
        tperiod = 2π / (abs(q2m) * 0.01)
        tspan = (0.0, 200 * tperiod)
        prob = ODEProblem(dummy_f, u0, tspan, param)

        init_Bmag = norm(gradient_B(u0[1:3], 0.0))
        dt_init = safety * 2π / (abs(q2m) * init_Bmag)

        sol = solve(prob, AdaptiveBoris(safety = safety); dt = dt_init)

        # Check if step sizes are not uniform (hence adapted)
        dts = diff(sol.t)
        @test !all(y -> isapprox(y, dts[1], rtol = 1.0e-5), dts)
        @test sol.t[end] ≈ tspan[2]

        # Energy conservation check
        total_energy(u) = 0.5 * m * norm(u[4:6])^2 + (q2m * m) * (-1.0e5 * u[2])
        E_start = total_energy(sol.u[1])
        E_end = total_energy(sol.u[end])
        @test isapprox(E_end, E_start, rtol = 1.0e-3)
    end

    @testset "Nonzero tspan[1]" begin
        # B field that is only present when t > 5
        B_field(r, t) = t > 5 ? SA[0.0, 0.0, 0.01] : SA[0.0, 0.0, 0.0]
        E_field(r, t) = SA[0.0, 0.0, 0.0]

        # param = (q2m, m, E, B, F)
        param = (-1.0e11, 1.0e-30, E_field, B_field, ZeroField())

        # Start at t=10. If absolute time is used, B should be 0.01.
        tspan = (10.0, 10.0 + 1.0e-7)
        dt = 1.0e-9
        u0 = SA[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
        prob = ODEProblem(dummy_f, u0, tspan, param)

        # Standard Boris
        sol_boris = solve(prob, Boris(); dt = dt)
        # If B was 0.01, vx should have deviated from initial 1.0e5
        @test abs(sol_boris.u[end][4] - 1.0e5) > 1.0e-4

        # Adaptive Boris
        sol_adaptive = solve(prob, AdaptiveBoris(safety = 0.1); dt = dt)
        @test abs(sol_adaptive.u[end][4] - 1.0e5) > 1.0e-4
    end

    @testset "SciML Output saving flags" begin
        # Setup
        x0 = [0.0, 0.0, 0.0]
        v0 = [0.0, 1.0e5, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 3.0e-8)
        dt = 3.0e-11

        zero_E(x, t) = SA[0.0, 0.0, 0.0]
        uniform_B(x, t) = SA[0.0, 0.0, 0.01]
        param = (-1.0e11, 1.0e-30, zero_E, uniform_B, zero_E)
        prob = ODEProblem(dummy_f, stateinit, tspan, param)

        # Baseline: save_everystep=true (default)
        sol = solve(prob, Boris(); dt)
        @test length(sol.t) == 1001 # 3e-8 / 3e-11 = 1000 steps + start

        # Scenario 2: Only final state
        sol = solve(prob, Boris(); dt, save_everystep = false, save_start = false, save_on = false)
        @test length(sol.t) == 1
        @test sol.t[1] ≈ tspan[2]

        # Scenario 3: Start and End
        sol = solve(prob, Boris(); dt, save_everystep = false, save_start = true, save_on = false)
        @test length(sol.t) == 2
        @test sol.t[1] == tspan[1]
        @test sol.t[end] ≈ tspan[2]

        # Multistep Boris test flags
        sol_ms = solve(prob, MultistepBoris2(n = 2); dt, save_everystep = false, save_start = true, save_on = false)
        @test length(sol_ms.t) == 2
        @test sol_ms.t[1] == tspan[1]
        @test sol_ms.t[end] ≈ tspan[2]
    end

    @testset "Custom parameter container" begin
        # The solvers must not depend on the layout of `p`, only on the accessors.
        param = CustomParam(1.0, 1.0, constant_Ey, constant_Bz)
        u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        prob = ODEProblem(dummy_f, u0, (0.0, 10.0), param)

        sol = solve(prob, Boris(); dt = 0.1)
        @test sol.u[end][1] ≈ 10.0 atol = 1.0e-6
        @test sol.u[end][4] ≈ 1.0 atol = 1.0e-6

        sol_ms = solve(prob, MultistepBoris4(n = 2); dt = 0.1)
        @test sol_ms.u[end][1] ≈ 10.0 atol = 1.0e-6
    end

    @testset "adaptive keyword" begin
        # `adaptive` is the switch between fixed and adaptive stepping, as it is
        # for any other SciML solver.
        param = (1.0, 1.0, ZeroField(), constant_Bz, ZeroField())
        u0 = SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        prob = ODEProblem(dummy_f, u0, (0.0, 2π), param)

        # Adaptive by default.
        sol = solve(prob, AdaptiveBoris(safety = 0.1); dt = 0.1)
        @test !all(y -> isapprox(y, diff(sol.t)[1], rtol = 1.0e-5), diff(sol.t))

        # An adaptive algorithm can be run with a fixed time step. The final step
        # is trimmed to land on the end of the time span, so it is excluded.
        sol = solve(prob, AdaptiveBoris(safety = 0.1); dt = 0.1, adaptive = false)
        dts = diff(sol.t)
        @test all(y -> isapprox(y, 0.1, rtol = 1.0e-10), dts[1:(end - 1)])
        @test sol.t[end] ≈ 2π

        # A fixed step algorithm refuses a request for adaptivity.
        @test_throws ArgumentError solve(prob, Boris(); dt = 0.1, adaptive = true)
    end
end
