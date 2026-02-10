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

    # 3. Dynamic GC ↔ FO Switching (Magnetic Bottle)
    #
    # B_z(z) = B0*(1 + α*z²), with B_r from ∇·B = 0:
    #   B_x ≈ -B0*α*x*z, B_y ≈ -B0*α*y*z
    #
    # Near the midplane (z ≈ 0) the field is weak and
    # curvature is high → FO. Away from midplane the
    # field strengthens and becomes more uniform → GC.
    # A bouncing particle triggers repeated GC → FO → GC
    # transitions.
    let
        B0 = 1.0e-4   # [T] background field
        α = 1.0e-4     # [m⁻²] mirror ratio parameter

        function bottle_B(x, t)
            Bz = B0 * (1 + α * x[3]^2)
            Bx = -B0 * α * x[1] * x[3]
            By = -B0 * α * x[2] * x[3]
            return SA[Bx, By, Bz]
        end
        B_bottle = TestParticle.Field(bottle_B)

        # Proton with moderate perpendicular velocity
        # and parallel velocity to bounce in the bottle
        x0 = SA[0.0, 0.0, 0.0]
        v_perp = 5.0e4   # [m/s]
        v_par = 2.0e5     # [m/s]
        v0 = SA[v_perp, 0.0, v_par]
        u0 = vcat(x0, v0)

        Ω = abs(q2m) * B0
        T_gyro = 2π / Ω
        tspan = (0.0, 50 * T_gyro)

        p = (
            q2m, m, E_field, B_bottle,
            TestParticle.ZeroField(),
        )
        alg = AdaptiveHybrid(;
            threshold = 0.05,
            dtmax = T_gyro,
            dtmin = 1.0e-4 * T_gyro,
        )

        sols = TestParticle.solve(TraceHybridProblem(u0, tspan, p), alg)
        sol = sols[1]
        @test sol.retcode == TestParticle.ReturnCode.Success
        @test length(sol.t) > 10

        # Energy conservation: kinetic energy should be
        # approximately conserved (no E field)
        KE_init = 0.5 * m * sum(abs2, v0)
        v_end = sol.u[end][SA[4, 5, 6]]
        KE_end = 0.5 * m * sum(abs2, v_end)
        @test KE_end ≈ KE_init rtol = 0.1

        # The particle should remain confined by the
        # mirror: z should not diverge
        z_vals = [u[3] for u in sol.u]
        @test maximum(abs, z_vals) < 1.0e6
    end
end
