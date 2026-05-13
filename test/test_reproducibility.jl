module test_reproducibility
using Test
using TestParticle
using OrdinaryDiffEq
using Random
using StaticArrays
using VelocityDistributionFunctions

@testset "Reproducibility" begin
    B_uniform(x) = SA[0, 0, 1.0e-8]
    E_zero = TestParticle.ZeroField()
    param = prepare(E_zero, B_uniform, species = Proton)
    tspan = (0.0, 1.0)
    stateinit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    prob = ODEProblem(trace!, stateinit, tspan, param)

    seed = 1234

    # Test OrdinaryDiffEq ensemble
    function prob_func(prob, ctx)
        v = rand(ctx.rng, SVector{3})
        remake(prob, u0 = [prob.u0[1:3]..., v...])
    end

    ensemble_prob = EnsembleProblem(prob; prob_func)

    sol1 = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories = 2, seed)
    sol2 = solve(ensemble_prob, Tsit5(), EnsembleSerial(); trajectories = 2, seed)

    @test sol1.u[1].u == sol2.u[1].u
    @test sol1.u[2].u == sol2.u[2].u

    # Test Boris ensemble
    trajs_boris1 = TestParticle.solve(TraceProblem(prob.u0, prob.tspan, prob.p), Boris(); dt = 0.1, trajectories = 2, seed)
    trajs_boris2 = TestParticle.solve(TraceProblem(prob.u0, prob.tspan, prob.p), Boris(); dt = 0.1, trajectories = 2, seed)

    @test trajs_boris1.u[1].u == trajs_boris2.u[1].u
    @test trajs_boris1.u[2].u == trajs_boris2.u[2].u

    # Test utility sampling
    rng1 = Xoshiro(seed)
    rng2 = Xoshiro(seed)
    @test sample_unit_sphere(rng1) == sample_unit_sphere(rng2)

    # sample_maxwellian is only available when VDF is loaded
    v1 = sample_maxwellian(rng1, 3000.0, TestParticle.mᵢ)
    v2 = sample_maxwellian(rng2, 3000.0, TestParticle.mᵢ)
    @test v1 == v2
end

end # module test_reproducibility
