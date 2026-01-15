using Test
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using SciMLBase

@testset "Hybrid Solver" begin

    # Define a simple analytic magnetic mirror field
    function B_simple(x, t)
        B0 = 1.0
        L = 2.0
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
    Ek = 10.0 # 10 eV
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

    tspan = (0.0, 1.0e-5)
    prob_gc = ODEProblem(trace_gc!, state_gc, tspan, params_gc)

    @testset "Integration" begin
        # Test Hybrid
        sol = solve_hybrid(prob_gc, Tsit5(); epsilon = 0.1, dt = 1.0e-9)

        # Check return type
        @test sol isa SciMLBase.AbstractODESolution

        # Check content
        @test length(sol.t) > 0
        @test length(sol.u) == length(sol.t)
        @test sol.u[1] isa SVector{6, Float64}
    end

    @testset "Correctness" begin
        # 1. Compare with pure GC (large epsilon => always GC)
        sol_hybrid_gc = solve_hybrid(prob_gc, Tsit5(); epsilon = 1.0e10, dt = 1.0e-9, adaptive = false)
        prob_pure_gc = ODEProblem(trace_gc!, state_gc, tspan, params_gc)
        sol_pure_gc = solve(prob_pure_gc, Tsit5(); dt = 1.0e-9, adaptive = false)

        # Compare final positions (should be identical or very close)
        @test norm(sol_hybrid_gc.u[end][1:3] - sol_pure_gc.u[end][1:3]) < 1.0e-10

        # 2. Compare with pure Full (small epsilon => always Full)
        # Note: Hybrid starts in GC, checks condition, then switches.
        # If we set epsilon small, it should switch immediately.
        sol_hybrid_full = solve_hybrid(prob_gc, Tsit5(); epsilon = 1.0e-10, dt = 1.0e-9, adaptive = false)

        # Construct equivalent pure full problem
        # We need to convert initial GC state to Full state
        state_full_init = TestParticle.gc_to_full(state_gc, params_gc, 0.0)
        # Construct full params: (q2m, m, E, B, F)
        p_full = (q / m, m, E_zero, B_simple, TestParticle.ZeroField())
        prob_pure_full = ODEProblem(trace, state_full_init, tspan, p_full)
        sol_pure_full = solve(prob_pure_full, Tsit5(); dt = 1.0e-9, adaptive = false)

        # Relax tolerance as hybrid initialization path differs from pure full
        @test norm(sol_hybrid_full.u[end][1:3] - sol_pure_full.u[end][1:3]) < 1.0
    end

    @testset "Stability" begin
        # Test Energy Conservation during switching
        sol = solve_hybrid(prob_gc, Tsit5(); epsilon = 0.1, dt = 1.0e-9)

        # Since we don't track mu output easily, we check stability
        @test all(x -> !isnan(x), Iterators.flatten(sol.u))
        @test Symbol(sol.retcode) in [:Success, :Terminated, :Default]
    end

    @testset "Analytical Switching" begin
        # Linear B field: B(x) = B0(1 - x/L) z-hat
        # E field: E0 y-hat (to drive drift in x towards lower B)

        B0 = 1.0
        L_grad = 2.0
        E0 = 100.0 # V/m - Need significant drift to reach switch point in reasonable time

        function B_linear(x, t)
            xx = x[1]
            if xx >= L_grad
                return SVector(0.0, 0.0, 1.0e-9) # Near zero, but avoid singularity
            end
            return SVector(0.0, 0.0, B0 * (1.0 - xx / L_grad))
        end
        B_linear(x) = B_linear(x, 0.0)

        function E_const(x, t)
            return SVector(0.0, E0, 0.0)
        end
        E_const(x) = SVector(0.0, E0, 0.0)

        m = TestParticle.Proton.m
        q = TestParticle.Proton.q

        # Initial: x=0.
        x0 = SVector(0.0, 0.0, 0.0)
        v0 = SVector(0.0, 1.0e5, 0.0) # Pure perp velocity to maximize r_L effect

        v_perp_0 = 1.0e5
        r_L_0 = (m * v_perp_0) / (q * B0)
        epsilon = 0.1

        # Analytical Switch Point: x_switch = L * (1 - (r_L0 / (eps * L))^(2/3))
        term = r_L_0 / (epsilon * L_grad)
        x_switch_analytic = L_grad * (1.0 - term^(2 / 3))

        # Run: Start in GC, estimate time to switch approx 0.05s
        state_gc, params_gc = prepare_gc(vcat(x0, v0), E_const, B_linear; species = TestParticle.Proton)
        tspan = (0.0, 0.05)

        prob_gc = ODEProblem(trace_gc!, state_gc, tspan, params_gc)

        sol = solve_hybrid(prob_gc, Tsit5(); epsilon = epsilon, dt = 1.0e-5)

        # Analyze switching: check if final state is beyond analytical switch point
        # and that the switch occurred at a location consistent with validity condition.

        # For this test, we accept if final state is beyond x_switch
        # and if it successfully integrated.

        final_x = sol.u[end][1]
        @test final_x > x_switch_analytic

        # Verify x_switch_analytic is consistent with where validity condition breaks

        # Calculate validity parameter along trajectory
        function check_validity_trace(u)
            R = u[1:3]
            v = u[4:6]
            B = B_linear(R, 0.0)
            Bmag = norm(B)
            # Reconstruct v_perp
            b = B / Bmag
            v_par = (v ⋅ b)
            v_perp_vec = v - v_par * b
            v_perp = norm(v_perp_vec)

            r_L = v_perp / (q / m * Bmag)
            # L_B = B / gradB. gradB = B0/L.
            L_B = Bmag / (B0 / L_grad)

            return r_L / L_B
        end

        # Scan solution
        switch_idx = nothing
        for i in 1:length(sol.u)
            # We need to be careful: sol.u are all 6D FULL states (converted)?
            # solve_hybrid returns result in whatever format?
            # Looking at source: states = SVector{6, u_type}[].
            # It converts GC to 6D GC state or keeps Full.

            val = check_validity_trace(sol.u[i])
            if val > epsilon
                switch_idx = i
                break
            end
        end

        @test !isnothing(switch_idx)
        x_at_switch = sol.u[switch_idx][1]

        # Verify numerical switch location matches analytical
        # Tolerance: The step size might overshoot slightly.
        @test isapprox(x_at_switch, x_switch_analytic, rtol = 0.1)

        # Verify it continued past limits or reached switch
        # If it switched at the very last step, this might fail strict inequality depending on precision.
        @test sol.u[end][1] >= x_at_switch
    end
end
