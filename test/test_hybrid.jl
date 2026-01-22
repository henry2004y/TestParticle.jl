using TestParticle
using StaticArrays
using LinearAlgebra
using Test

@testset "Adaptive Hybrid Solver" begin
    # 1. Pure Adiabatic Case (Uniform Field)
    B_func(x, t) = SA[0.0, 0.0, 0.01]
    B_field = TestParticle.Field(B_func)
    E_field = TestParticle.Field((x, t) -> SA[0.0, 0.0, 0.0])

    m = TestParticle.mᵢ
    q = TestParticle.qᵢ
    q2m = q / m

    x0 = SA[0.0, 0.0, 0.0]
    v0 = SA[1.0e4, 0.0, 1.0e3]
    u0 = vcat(x0, v0)
    tspan = (0.0, 1.0e-4)

    p = (q2m, m, E_field, B_field)
    alg = AdaptiveHybrid(threshold = 0.1, dtmax = 1.0e-6)

    sols = TestParticle.solve(TraceHybridProblem(u0, tspan, p), alg)
    @test sols[1].retcode == TestParticle.ReturnCode.Success

    # 2. Switching Case (Curved/Sheared Field)
    # B = [B0*cos(kx), B0*sin(kx), 0]
    # This field line is curved in the xy plane.
    function sheared_B_func(x, t)
        B0 = 0.01
        k = 100.0 # High curvature
        return SA[B0 * cos(k * x[1]), B0 * sin(k * x[1]), 0.0]
    end
    sheared_B = TestParticle.Field(sheared_B_func)

    p2 = (q2m, m, E_field, sheared_B)
    # Threshold is 0.1. ρ ≈ 0.01. Rc ≈ 1/k = 0.01.
    # ϵ = ρ / Rc ≈ 1.0 > 0.1 => Should start in FO or switch to FO.

    tspan2 = (0.0, 1.0e-5)
    sols2 = TestParticle.solve(TraceHybridProblem(u0, tspan2, p2), alg)
    @test sols2[1].retcode == TestParticle.ReturnCode.Success

    # Let's use a lower k for adiabatic start and higher k for transition if possible.
    # Or just verify with logs/print.
end
