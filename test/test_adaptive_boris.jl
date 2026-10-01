module test_adaptive_boris
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
        @test alg isa AdaptiveBoris
        @test alg.safety == 0.2

        alg_def = AdaptiveBoris()
        @test alg_def.safety == 0.1

        # Fixed and adaptive stepping are distinct types, and `adaptive` is the
        # keyword that selects between them at solve time.
        alg_plain = Boris()
        @test alg_plain isa Boris
        @test !(alg_plain isa AdaptiveBoris)
    end

    @testset "Solve Integration" begin
        @testset "Single trajectory" begin
            sol = TestParticle.solve(
                prob, AdaptiveBoris()
            )
            @test sol.retcode ==
                TestParticle.ReturnCode.Success
            @test length(sol.t) > 1
        end

        @testset "EnsembleProblem" begin
            # A TraceProblem goes into a SciML ensemble as it is; SciML calls
            # `solve` on each trajectory it builds.
            trajectories = 4
            sols = TestParticle.solve(
                EnsembleProblem(prob), AdaptiveBoris(), EnsembleThreads();
                trajectories
            )
            @test length(sols.u) == trajectories
            @test all(s -> s.retcode == TestParticle.ReturnCode.Success, sols.u)
        end

        @testset "EnsembleThreads" begin
            trajectories = 4
            sols = TestParticle.solve(
                prob, AdaptiveBoris(), EnsembleThreads();
                trajectories = trajectories
            )
            @test length(sols.u) == trajectories
            @test all(s -> s.retcode == TestParticle.ReturnCode.Success, sols.u)
        end
    end
end

end # module test_adaptive_boris
