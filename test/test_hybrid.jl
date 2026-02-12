using TestParticle
using StaticArrays
using LinearAlgebra
using Test

@testset "Adaptive Hybrid Solver" begin
    m = TestParticle.mᵢ
    q = TestParticle.qᵢ
    q2m = q / m

    E_field = TestParticle.Field((x, t) -> SA[0.0, 0.0, 0.0])

    # 1. Pure Adiabatic Case (Uniform Field)
    let
        B_func(x, t) = SA[0.0, 0.0, 0.01]
        B_field = TestParticle.Field(B_func)

        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[1.0e4, 0.0, 1.0e3]
        u0 = vcat(x0, v0)
        tspan = (0.0, 1.0e-4)

        p = (q2m, m, E_field, B_field, TestParticle.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)

        sols = TestParticle.solve(TraceHybridProblem(u0, tspan, p), alg)
        @test sols[1].retcode == TestParticle.ReturnCode.Success
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
        sheared_B = TestParticle.Field(sheared_B_func)

        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[1.0e4, 0.0, 1.0e3]
        u0 = vcat(x0, v0)
        tspan = (0.0, 1.0e-5)

        p = (q2m, m, E_field, sheared_B, TestParticle.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)

        sols = TestParticle.solve(TraceHybridProblem(u0, tspan, p), alg)
        @test sols[1].retcode == TestParticle.ReturnCode.Success
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
        B_bottle = TestParticle.Field(bottle_B)

        x0 = SA[0.0, 0.0, 0.0]
        v_perp = 5.0e4  # [m/s]
        v_par = 1.0e5   # [m/s]
        v0 = SA[v_perp, 0.0, v_par]
        u0 = vcat(x0, v0)

        Ω = abs(q2m) * B0
        T_gyro = 2π / Ω
        tspan = (0.0, 30 * T_gyro)

        p = (q2m, m, E_field, B_bottle, TestParticle.ZeroField())

        # Consistent with demo_hybrid.jl
        alg = AdaptiveHybrid(;
            threshold = 0.1,
            dtmax = T_gyro,
            dtmin = 1.0e-4 * T_gyro,
            check_interval = 100,
        )

        sols = TestParticle.solve(TraceHybridProblem(u0, tspan, p), alg)

        @test sols[1].retcode == TestParticle.ReturnCode.Success
        # Check that we actually have a reasonable number of steps (hybrid should adapt)
        @test length(sols[1].t) > 100

        # Energy conservation check (approximate for hybrid)
        KE_init = 0.5 * m * sum(abs2, v0)
        KE_final = 0.5 * m * sum(abs2, sols[1].u[end][SA[4, 5, 6]])
        @test isapprox(KE_final, KE_init; rtol = 0.1)
    end
end
