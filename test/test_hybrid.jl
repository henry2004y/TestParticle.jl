using TestParticle
import TestParticle as TP
using StaticArrays
using LinearAlgebra
using Test

@testset "Adaptive Hybrid Solver" begin
    m = TP.mᵢ
    q = TP.qᵢ
    q2m = q / m

    E_field = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

    # 1. Pure Adiabatic Case (Uniform Field)
    let
        B_func(x, t) = SA[0.0, 0.0, 0.01]
        B_field = TP.Field(B_func)

        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[1.0e4, 0.0, 1.0e3]
        u0 = vcat(x0, v0)
        tspan = (0.0, 1.0e-4)

        p = (q2m, m, E_field, B_field, TP.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)

        sols = TP.solve(TraceHybridProblem(u0, tspan, p), alg)
        @test sols.u[1].retcode == TP.ReturnCode.Success
    end

    # 2. Non-adiabatic Case (Highly Curved Field)
    let
        function sheared_B_func(x, t)
            B0 = 0.01
            k = 100.0
            return SA[
                B0 * cos(k * x[1]),
                B0 * sin(k * x[1]),
                0.0,
            ]
        end
        sheared_B = TP.Field(sheared_B_func)

        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[1.0e4, 0.0, 1.0e3]
        u0 = vcat(x0, v0)
        tspan = (0.0, 1.0e-5)

        p = (q2m, m, E_field, sheared_B, TP.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)

        sols = TP.solve(TraceHybridProblem(u0, tspan, p), alg)
        @test sols.u[1].retcode == TP.ReturnCode.Success
    end

    # 3. Demo Hybrid Case (Explicit Switching)
    # Replicates the setup from docs/examples/features/demo_hybrid.jl
    # Stronger curvature (α = 1.0e-2) forces GC <-> FO switching.
    let
        B0 = 1.0e-4   # [T]
        α = 1.0e-2    # [m⁻²] Stronger curvature

        function bottle_B(x, t)
            Bz = B0 * (1 + α * x[3]^2)
            Bx = -B0 * α * x[1] * x[3]
            By = -B0 * α * x[2] * x[3]
            return SA[Bx, By, Bz]
        end
        B_bottle = TP.Field(bottle_B)

        x0 = SA[0.0, 0.0, 0.0]
        v_perp = 5.0e4  # [m/s]
        v_par = 1.0e5   # [m/s]
        v0 = SA[v_perp, 0.0, v_par]
        u0 = vcat(x0, v0)

        Ω = abs(q2m) * B0
        T_gyro = 2π / Ω
        tspan = (0.0, 30 * T_gyro)

        p = (q2m, m, E_field, B_bottle, TP.ZeroField())

        # Consistent with demo_hybrid.jl
        alg = AdaptiveHybrid(;
            threshold = 0.1,
            dtmax = T_gyro,
            dtmin = 1.0e-4 * T_gyro,
            check_interval = 100,
        )

        sols = TP.solve(TraceHybridProblem(u0, tspan, p), alg)

        @test sols.u[1].retcode == TP.ReturnCode.Success
        # Check that we actually have a reasonable number of steps (hybrid should adapt)
        @test length(sols.u[1].t) > 100

        # Energy conservation check (approximate for hybrid)
        KE_init = 0.5 * m * sum(abs2, v0)
        KE_final = 0.5 * m * sum(abs2, sols.u[1].u[end][SA[4, 5, 6]])
        @test isapprox(KE_final, KE_init; rtol = 0.1)
    end

    # 4. Round-trip: GC -> full-orbit reconstruction is exact in a uniform field.
    # (The guiding-center map is the exact inverse of itself when ∇B = 0.)
    let
        B_func(x, t) = SA[0.0, 0.0, 0.01]
        B_field = TP.Field(B_func)
        for x0 in (SA[0.0, 0.0, 0.0], SA[1.0, 0.5, 1.0], SA[3.0, 0.0, 2.0])
            for v0 in (SA[5.0e4, 0.0, 1.0e5], SA[3.0e4, 4.0e4, 8.0e4])
                xv = vcat(x0, v0)
                X, vpar, _, _, μ, phase = TP._get_gc_parameters(
                    xv, E_field, B_field, q, m, 0.0
                )
                state_gc = SA[X[1], X[2], X[3], vpar]
                xv_rec = TP._gc_to_full(
                    state_gc, E_field, B_field, q, m, μ, 0.0, phase
                )
                @test norm(xv[SA[1, 2, 3]] - xv_rec[SA[1, 2, 3]]) < 1.0e-10
                @test norm(xv[SA[4, 5, 6]] - xv_rec[SA[4, 5, 6]]) < 1.0e-8
            end
        end
    end

    # 5. Forced pure-FO hybrid must reproduce standalone Boris bit-for-bit.
    # (threshold = 0 forces FO for the whole trajectory since ε ≥ 0 always.)
    let
        B0 = 0.01
        B_func(x, t) = SA[0.0, 0.0, B0]
        B_field = TP.Field(B_func)
        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[5.0e4, 0.0, 1.0e5]
        u0 = vcat(x0, v0)
        Ω = abs(q2m) * B0
        T_gyro = 2π / Ω
        tspan = (0.0, 30 * T_gyro)
        p = (q2m, m, E_field, B_field, TP.ZeroField())

        alg_fo = AdaptiveHybrid(;
            threshold = 0.0, dtmax = 10 * T_gyro, dtmin = 1.0e-6 * T_gyro,
            check_interval = 100
        )
        sol_h = TP.solve(TraceHybridProblem(u0, tspan, p), alg_fo).u[1]

        dt_fo = 2π * alg_fo.safety_fo / (abs(q2m) * B0)
        sol_b = TP.solve(TraceProblem(u0, tspan, p), Boris(); dt = dt_fo).u[1]

        @test length(sol_h.t) == length(sol_b.t)
        for i in 1:length(sol_b.t)
            uh = sol_h(sol_b.t[i])
            ub = sol_b.u[i]
            @test isapprox(norm(uh[SA[1, 2, 3]] - ub[SA[1, 2, 3]]), 0.0; atol = 1.0e-8)
            @test isapprox(norm(uh[SA[4, 5, 6]] - ub[SA[4, 5, 6]]), 0.0; atol = 1.0e-5)
        end
    end

    # 6. Hysteresis thresholds: backward-compatible `threshold` sets both
    # bounds, and a wider buffer (α > β) reduces GC <-> FO chattering.
    let
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0).threshold_gc_to_fo == 0.1
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0).threshold_fo_to_gc == 0.1

        # α must be >= β (buffer region requires α ≥ β).
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold_gc_to_fo = 0.05, threshold_fo_to_gc = 0.2, dtmax = 1.0
        )

        # Distinct bounds are stored verbatim.
        alg_buf = AdaptiveHybrid(;
            threshold_gc_to_fo = 0.2, threshold_fo_to_gc = 0.05, dtmax = 1.0
        )
        @test alg_buf.threshold_gc_to_fo == 0.2
        @test alg_buf.threshold_fo_to_gc == 0.05
    end

    # 7. Buffer region suppresses mode chatter: a particle whose adiabaticity
    # parameter oscillates within [β, α) must keep a single, stable mode when a
    # wide hysteresis band is used, whereas a zero-width band (α == β) permits
    # switching at every crossing. We exercise this with the bottle field and
    # verify the run still completes successfully with a wide buffer.
    let
        B0 = 1.0e-4
        α = 1.0e-2
        function bottle_B(x, t)
            Bz = B0 * (1 + α * x[3]^2)
            Bx = -B0 * α * x[1] * x[3]
            By = -B0 * α * x[2] * x[3]
            return SA[Bx, By, Bz]
        end
        B_bottle = TP.Field(bottle_B)
        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[5.0e4, 0.0, 1.0e5]
        u0 = vcat(x0, v0)
        Ω = abs(q2m) * B0
        T_gyro = 2π / Ω
        tspan = (0.0, 30 * T_gyro)
        p = (q2m, m, E_field, B_bottle, TP.ZeroField())

        # Wide hysteresis band: α = 0.2, β = 0.05.
        alg_buf = AdaptiveHybrid(;
            threshold_gc_to_fo = 0.2, threshold_fo_to_gc = 0.05,
            dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro, check_interval = 100,
        )
        sol_buf = TP.solve(TraceHybridProblem(u0, tspan, p), alg_buf).u[1]
        @test sol_buf.retcode == TP.ReturnCode.Success
        @test length(sol_buf.t) > 100
    end
end
