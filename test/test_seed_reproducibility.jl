module test_seed_reproducibility

using Test
using TestParticle
using StaticArrays
using OrdinaryDiffEq
using Random
using LinearAlgebra
using VelocityDistributionFunctions

# Set up simple fields that support ForwardDiff seamlessly
B_uniform(x) = SA[0.0, 0.0, 1.0e-8]
E_uniform(x) = SA[0.0, 0.0, 0.0]
tspan = (0.0, 1.0e-4)

# Create guiding center parameters and initial state
x0 = SA[0.0, 0.0, 0.0]
v0 = SA[1.0e4, 0.0, 1.0e4]
u0_full = vcat(x0, v0)
u0_gc, param = TestParticle.prepare_gc(
    u0_full, E_uniform, B_uniform; species = Proton
)

function prob_func(prob, ctx)
    vdf = TestParticle.Maxwellian(SA[0.0, 0.0, 0.0], 1.0)
    # Using ctx.rng is critical here!
    v = rand(ctx.rng, vdf)
    # Rk4/Rk45 expects a 4-element state (x, y, z, vpar)
    return remake(prob; u0 = [prob.u0[1:3]..., v[3]])
end

@testset "GC Solver Seed Reproducibility" begin
    let
        prob = TraceGCProblem(u0_gc, tspan, param; prob_func)
        seed = 1234

        # Test RK4
        sols1_rk4 = TestParticle.solve(
            prob; trajectories = 3, dt = 1.0e-5, alg = :rk4, seed
        )
        sols2_rk4 = TestParticle.solve(
            prob; trajectories = 3, dt = 1.0e-5, alg = :rk4, seed
        )

        @test length(sols1_rk4.u) == 3
        for i in 1:3
            @test sols1_rk4.u[i].u == sols2_rk4.u[i].u
        end

        # Test RK45
        sols1_rk45 = TestParticle.solve(
            prob; trajectories = 3, dt = 1.0e-5, alg = :rk45, seed
        )
        sols2_rk45 = TestParticle.solve(
            prob; trajectories = 3, dt = 1.0e-5, alg = :rk45, seed
        )

        @test length(sols1_rk45.u) == 3
        for i in 1:3
            @test sols1_rk45.u[i].u == sols2_rk45.u[i].u
        end
    end
end

@testset "Hybrid Solver Seed Reproducibility" begin
    let
        m = TestParticle.mᵢ
        q = TestParticle.qᵢ
        q2m = q / m
        E_zero = TestParticle.Field((x, t) -> SA[0.0, 0.0, 0.0])

        # Highly curved field to force switching which calls rand(rng)
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

        p = (q2m, m, E_zero, sheared_B, TestParticle.ZeroField())

        # Setup prob_func that samples initial conditions randomly
        function hybrid_prob_func(prob, ctx)
            # Use rand(ctx.rng) to randomize velocity slightly
            v = SA[1.0e4 + 1.0e2 * rand(ctx.rng), 0.0, 1.0e3]
            # Ensure u0 is constructed as a static SVector{6, Float64}
            u0_new = SVector{6, Float64}(
                prob.u0[1], prob.u0[2], prob.u0[3],
                v[1], v[2], v[3]
            )
            return remake(prob; u0 = u0_new)
        end

        prob = TraceHybridProblem(u0, tspan, p; prob_func = hybrid_prob_func)
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)
        seed = 5678

        sols1 = TestParticle.solve(prob, alg; trajectories = 3, seed)
        sols2 = TestParticle.solve(prob, alg; trajectories = 3, seed)

        @test length(sols1.u) == 3
        for i in 1:3
            @test sols1.u[i].u == sols2.u[i].u
        end
    end
end

@testset "Seed Independence (Boris / GC / Hybrid)" begin
    # Distinct seeds must give *independent* ensembles, not a rotation of the same
    # RNG stream. (A naive `Xoshiro(seed + i)` only shifts which RNG each trajectory
    # index draws from, so seed=2 was previously a near-rotation of seed=1.)

    # --- Boris (full 6D) ---
    let
        param_b = prepare(E_uniform, B_uniform; species = Proton)
        function pf_b(prob, ctx)
            v = SA[1.0e4 * rand(ctx.rng), 1.0e4 * rand(ctx.rng), 1.0e4 * rand(ctx.rng)]
            return remake(prob; u0 = SA[0.0, 0.0, 0.0, v...])
        end
        prob = TraceProblem(SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param_b; prob_func = pf_b)
        init_v(sols) = [
            SVector(sols.u[i].u[1][4], sols.u[i].u[1][5], sols.u[i].u[1][6])
                for i in 1:length(sols.u)
        ]
        v1 = init_v(TestParticle.solve(prob, Boris(); trajectories = 6, dt = 1.0e-6, seed = 1))
        v2 = init_v(TestParticle.solve(prob, Boris(); trajectories = 6, dt = 1.0e-6, seed = 2))
        v1b = init_v(TestParticle.solve(prob, Boris(); trajectories = 6, dt = 1.0e-6, seed = 1))
        @test v1 == v1b
        @test !all(v2[i] == v1[i + 1] for i in 1:5)
    end

    # --- GC (rk4) ---
    let
        prob = TraceGCProblem(u0_gc, tspan, param; prob_func)
        init_v(sols) = [sols.u[i].u[1][4] for i in 1:length(sols.u)]
        v1 = init_v(TestParticle.solve(prob; trajectories = 6, dt = 1.0e-5, alg = :rk4, seed = 1))
        v2 = init_v(TestParticle.solve(prob; trajectories = 6, dt = 1.0e-5, alg = :rk4, seed = 2))
        v1b = init_v(TestParticle.solve(prob; trajectories = 6, dt = 1.0e-5, alg = :rk4, seed = 1))
        @test v1 == v1b
        @test !all(v2[i] == v1[i + 1] for i in 1:5)
    end

    # --- Hybrid ---
    let
        m = TestParticle.mᵢ
        q = TestParticle.qᵢ
        q2m = q / m
        E_zero = TestParticle.ZeroField()
        function sheared_B_func(x, t)
            B0 = 0.01
            k = 100.0
            return SA[B0 * cos(k * x[1]), B0 * sin(k * x[1]), 0.0]
        end
        sheared_B = TestParticle.Field(sheared_B_func)
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        p = (q2m, m, E_zero, sheared_B, TestParticle.ZeroField())
        function pf_h(prob, ctx)
            v = SA[1.0e4 + 1.0e2 * rand(ctx.rng), 0.0, 1.0e3]
            return remake(
                prob; u0 = SVector{6, Float64}(
                    prob.u0[1], prob.u0[2], prob.u0[3], v[1], v[2], v[3]
                )
            )
        end
        prob = TraceHybridProblem(u0, tspan, p; prob_func = pf_h)
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)
        init_v(sols) = [sols.u[i].u[1][4] for i in 1:length(sols.u)]
        v1 = init_v(TestParticle.solve(prob, alg; trajectories = 6, seed = 1))
        v2 = init_v(TestParticle.solve(prob, alg; trajectories = 6, seed = 2))
        v1b = init_v(TestParticle.solve(prob, alg; trajectories = 6, seed = 1))
        @test v1 == v1b
        @test !all(v2[i] == v1[i + 1] for i in 1:5)
    end
end

end # module test_seed_reproducibility
