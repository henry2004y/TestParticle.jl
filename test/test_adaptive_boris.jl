using TestParticle
using StaticArrays
using LinearAlgebra
using Statistics
using Test

@testset "Adaptive Boris" begin
    # 1. Test basic functionality with gradient field (variable step)
    function gradient_B_field(x)
        # B = (0, 0, 1 + z^2)
        return SVector{3}(0.0, 0.0, 1.0 + x[3]^2)
    end

    function zero_E_field(x)
        return SVector{3}(0.0, 0.0, 0.0)
    end

    u0 = [1.0, 0.0, 0.0, 0.1, 0.1, 1.0] # x, y, z, vx, vy, vz
    tspan = (0.0, 10.0)
    # Use prepare to generate compatible parameter tuple
    # prepare(E, B; species=Proton) -> (q2m, m, E, B, F)
    # We want q2m=1. We can set q=1, m=1.
    param = prepare(zero_E_field, gradient_B_field; species=User, q=1.0, m=1.0)
    prob = TraceProblem(u0, tspan, param)

    # Run fixed step for comparison
    sol_fixed = TestParticle.solve(prob; dt=0.1)[1]

    # Run adaptive step
    alg = AdaptiveBoris(dtmax=0.5, dtmin=1e-4, safety=0.1)
    sol_adaptive = TestParticle.solve(prob, alg; dt=0.1)[1]

    steps = diff(sol_adaptive.t)

    @test length(sol_adaptive.t) > length(sol_fixed.t) # Expect more steps as B increases and dt decreases
    @test maximum(steps) <= 0.5
    # @test minimum(steps) >= 1e-4
    # The last step might be very small to hit tspan[2] exactly.
    # We should check if steps EXCEPT the last one obey the limit, or just relax this check.
    # If dt_next is clamped, it should be >= dtmin.
    # However, for the FINAL step, we do `dt_next = tspan[2] - t_current`, which can be arbitrarily small.
    # So we should filter out the last step or check `minimum(steps[1:end-1])`.

    if length(steps) > 1
        @test minimum(steps[1:end-1]) >= 1e-4
    end
    @test !all(isapprox.(steps, steps[1])) # Steps should vary

    # Check that it saves correctly
    @test sol_adaptive.t[1] == 0.0
    @test sol_adaptive.t[end] == 10.0

    # 2. Test with constant field (should reach dtmax)
    function constant_B_field(x)
        return SVector{3}(0.0, 0.0, 1.0)
    end

    param_const = prepare(zero_E_field, constant_B_field; species=User, q=1.0, m=1.0)
    prob_const = TraceProblem(u0, tspan, param_const)

    # q2m=1, B=1. safety=0.1 => dt_ideal = 0.1.
    # Set dtmax=0.05. It should clamp to 0.05.
    alg_const = AdaptiveBoris(dtmax=0.05, dtmin=1e-4, safety=0.1)
    sol_const = TestParticle.solve(prob_const, alg_const; dt=0.01)[1]

    steps_const = diff(sol_const.t)
    # The first step uses initial dt=0.01. Subsequent steps should adapt.
    # Since dt_ideal=0.1 and dtmax=0.05, it should use 0.05.
    # Wait, the first step in the loop calculates dt based on B.
    # The first `update_velocity!` push back uses initial dt=0.01.
    # Then loop starts. `dt_next` is calculated.
    # So `t[2] - t[1]` should be `dt_next`.

    # @test steps_const[2] ≈ 0.05 atol=1e-4
    # Note: steps_const[1] is t[2]-t[1].
    # In my implementation:
    # `dt_next` is calculated at start of loop.
    # `update_velocity!` uses `dt_next`.
    # `push!` saves `t_current + dt_next`.
    # So t[2] = t[1] + dt_next.
    # Thus steps_const[1] should be dt_next (approx 0.05).

    @test isapprox(steps_const[2], 0.05, atol=1e-4)
    @test all(s -> s <= 0.05 + 1e-9, steps_const)

    # 3. Test exact end time
    @test sol_const.t[end] == 10.0
end
