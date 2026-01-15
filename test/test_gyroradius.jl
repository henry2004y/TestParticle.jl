using Test
using StaticArrays
using TestParticle
using OrdinaryDiffEq

@testset "Gyroradius Utility" begin
    @testset "Non-relativistic" begin
        # Uniform B field
        B0 = 1.0e-8
        B_func(x, t) = SA[0.0, 0.0, B0]
        E_func(x, t) = SA[0.0, 0.0, 0.0]

        # Particle parameters
        q = TestParticle.qᵢ
        m = TestParticle.mᵢ

        # Initial state: v perpendicular to B
        v_perp = 1.0e5
        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[v_perp, 0.0, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 1.0)

        param = prepare(E_func, B_func; species = Ion(1, 1)) # Proton
        prob = ODEProblem(trace!, stateinit, tspan, param)
        sol = solve(prob, Tsit5())

        # Theoretical gyroradius
        # r = m * v_perp / (q * B)
        r_expected = m * v_perp / (q * B0)

        r_calc = get_gyroradius(sol, 0.5)
        @test r_calc ≈ r_expected rtol = 1.0e-5
    end

    @testset "Relativistic" begin
        # Uniform B field
        B0 = 0.01
        B_func(x, t) = SA[0.0, 0.0, B0]
        E_func(x, t) = SA[0.0, 0.0, 0.0]

        # Particle parameters
        q = TestParticle.qₑ
        m = TestParticle.mₑ
        c = TestParticle.c

        # Initial state: relativistic velocity
        # Let v = 0.5c
        # γ = 1/sqrt(1-0.5^2) = 1/sqrt(0.75) ≈ 1.1547
        # u = γv
        v_mag = 0.5 * c
        gamma = 1.0 / sqrt(1 - (v_mag / c)^2)
        u_perp = gamma * v_mag

        x0 = SA[0.0, 0.0, 0.0]
        u0 = SA[u_perp, 0.0, 0.0] # γv
        stateinit = [x0..., u0...]
        tspan = (0.0, 1.0e-9)

        param = prepare(E_func, B_func; species = Electron)
        prob = ODEProblem(trace_relativistic!, stateinit, tspan, param)
        sol = solve(prob, Vern6())

        # Theoretical relativistic gyroradius
        # r = p_perp / (q * B) = m * (γv)_perp / (q * B)
        r_expected = m * u_perp / (abs(q) * B0)

        r_calc = get_gyroradius(sol, 0.5e-9)
        @test r_calc ≈ r_expected rtol = 1.0e-3
    end

    @testset "Zero Field" begin
        # B = 0 -> r = Inf
        B_func(x, t) = SA[0.0, 0.0, 0.0]
        E_func(x, t) = SA[0.0, 0.0, 0.0]

        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[1.0, 0.0, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 1.0)

        param = prepare(E_func, B_func; species = Ion(1, 1))
        prob = ODEProblem(trace!, stateinit, tspan, param)
        sol = solve(prob, Tsit5())

        @test get_gyroradius(sol, 0.5) == Inf
    end

    @testset "ExB Drift Correction" begin
        # B = 10 nT, E in y direction
        B0 = 10.0e-9
        E0 = 1.0e-4 # V/m
        B_func(x, t) = SA[0.0, 0.0, B0]
        E_func(x, t) = SA[0.0, E0, 0.0]

        # Drift velocity: v_E = (E x B) / B^2
        # E = (0, E0, 0), B = (0, 0, B0)
        # E x B = (E0*B0, 0, 0)
        # v_E = (E0/B0, 0, 0)
        v_drift = E0 / B0 # 1e-4 / 10e-9 = 10,000 m/s

        # Initial state: particle with exactly the drift velocity
        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[v_drift, 0.0, 0.0]
        stateinit = [x0..., v0...]
        tspan = (0.0, 1.0)

        param = prepare(E_func, B_func; species = Ion(1, 1))
        prob = ODEProblem(trace!, stateinit, tspan, param)
        sol = solve(prob, Tsit5())

        # With pure drift, gyroradius should be effectively 0
        r_calc = get_gyroradius(sol, 0.5)
        @test r_calc ≈ 0.0 atol = 1.0e-3
    end

    @testset "Vector Dispatch" begin
        q = TestParticle.qᵢ
        m = TestParticle.mᵢ
        B0 = 1.0e-8
        v_perp = 1.0e5

        # Case 1: Vector inputs for V and B, same result as scalar
        V = SA[v_perp, 0.0, 0.0]
        B = SA[0.0, 0.0, B0]
        r_scalar = get_gyroradius(v_perp, B0; q, m)
        r_vector = get_gyroradius(V, B; q, m)
        @test r_vector ≈ r_scalar

        # Case 2: Vector inputs with E field correction
        E0 = 1.0e-4
        E = SA[0.0, E0, 0.0]
        B_vec = SA[0.0, 0.0, 10.0e-9] # same as ExB drift correction test
        v_drift = E0 / 10.0e-9
        V_drift_vec = SA[v_drift, 0.0, 0.0]

        # Pure drift should give ~0 gyroradius
        r_drift = get_gyroradius(V_drift_vec, B_vec, E; q, m)
        @test r_drift ≈ 0.0 atol = 1.0e-3

        # Compare with explicit scalar calculation without E field (which would be wrong)
        # Without E, r would be m*v_drift / qB
        r_wrong = get_gyroradius(v_drift, 10.0e-9; q, m)
        @test r_wrong > 0.1 # Should be large
    end
end
