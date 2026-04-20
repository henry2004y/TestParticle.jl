module test_adaptive_multistep_boris

using Test
using TestParticle
using StaticArrays

@testset "Adaptive Multistep Boris" begin
    B_field(r, t) = SA[0, 0, 1.0e-8]
    E_field(r, t) = SA[0, 0, 0]
    param = TestParticle.prepare(E_field, B_field)

    tspan = (0.0, 1.0)
    x0 = [1.0, 0.0, 0.0]
    v0 = [0.0, 1.0, 0.0]
    stateinit = [x0..., v0...]
    prob = TestParticle.TraceProblem(stateinit, tspan, param)

    @testset "Constructor" begin
        # Test AdaptiveMultistepBoris constructors
        alg2 = AdaptiveMultistepBoris{2}(n = 1, safety = 0.2)
        @test alg2 isa MultistepBoris{2, true}
        @test alg2.n == 1
        @test alg2.safety == 0.2

        alg4 = AdaptiveMultistepBoris{4}(n = 2)
        @test alg4 isa MultistepBoris{4, true}
        @test alg4.n == 2
        @test alg4.safety == 0.1

        alg6 = AdaptiveMultistepBoris{6}()
        @test alg6 isa MultistepBoris{6, true}
        @test alg6.n == 1
        @test alg6.safety == 0.1

        # Test MultistepBoris default still works (fixed-step)
        alg_fixed = MultistepBoris{2}(n = 1)
        @test alg_fixed isa MultistepBoris{2, false}
        @test alg_fixed.safety == 0.0
    end

    @testset "Consistency with AdaptiveBoris" begin
        # AdaptiveMultistepBoris{2}(n=1) should be identical to AdaptiveBoris
        sol_std = TestParticle.solve(prob, AdaptiveBoris(safety = 0.1))[1]
        sol_multi = TestParticle.solve(prob, AdaptiveMultistepBoris{2}(n = 1, safety = 0.1))[1]

        @test sol_std.t ≈ sol_multi.t
        @test sol_std.u ≈ sol_multi.u
    end

    @testset "Higher Order Integration" begin
        # Test N=4 and N=6
        trajectories = 2
        for N in [4, 6]
            sols = TestParticle.solve(
                prob, AdaptiveMultistepBoris{N}(n = 2), EnsembleSerial(); trajectories
            )
            @test length(sols) == trajectories
            @test all(s -> s.retcode == TestParticle.ReturnCode.Success, sols)
        end
    end

    @testset "Ensemble Algorithms" begin
        trajectories = 4
        alg = AdaptiveMultistepBoris{4}(n = 1)

        @testset "EnsembleThreads" begin
            sols = TestParticle.solve(
                prob, alg, EnsembleThreads(); trajectories
            )
            @test length(sols) == trajectories
            @test all(s -> s.retcode == TestParticle.ReturnCode.Success, sols)
        end

        @testset "EnsembleDistributed" begin
            # Note: EnsembleDistributed might be slow in tests, but we should at least check if it dispatches
            sols = TestParticle.solve(
                prob, alg, EnsembleDistributed(); trajectories
            )
            @test length(sols) == trajectories
            @test all(s -> s.retcode == TestParticle.ReturnCode.Success, sols)
        end
    end
end

end # module test_adaptive_multistep_boris
