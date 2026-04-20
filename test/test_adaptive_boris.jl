using Test
using TestParticle
using StaticArrays

@testset "AdaptiveBoris Compatibility" begin
    B_field(r, t) = SA[0, 0, 1.0e-8]
    E_field(r, t) = SA[0, 0, 0]
    param = TestParticle.prepare(E_field, B_field)

    tspan = (0.0, 1.0)
    x0 = [1.0, 0.0, 0.0]
    v0 = [0.0, 1.0, 0.0]
    stateinit = [x0..., v0...]
    prob = TestParticle.TraceProblem(
        stateinit, tspan, param
    )

    @testset "Constructor" begin
        alg = AdaptiveBoris(; safety = 0.2)
        @test alg isa Boris{true}
        @test alg isa AdaptiveBoris
        @test alg.safety == 0.2

        alg_def = AdaptiveBoris()
        @test alg_def.safety == 0.1

        alg_plain = Boris()
        @test alg_plain isa Boris{false}
        @test alg_plain.safety == 0.0
    end

    @testset "Solve Integration" begin
        sols = TestParticle.solve(
            prob, AdaptiveBoris()
        )
        @test sols[1].retcode ==
            TestParticle.ReturnCode.Success
        @test length(sols[1].t) > 1
    end
end
