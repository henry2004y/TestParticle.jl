module test_boris_gpu

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using KernelAbstractions
using LinearAlgebra

uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
zero_E(x) = SA[0.0, 0.0, 0.0]

"The largest difference between two sets of states, relative to the state itself."
function max_rel_diff(a, b)
    return maximum(zip(a, b)) do (u, v)
        return norm(collect(u) - collect(v)) / max(norm(collect(v)), one(eltype(v)))
    end
end

@testset "Boris on a backend" begin
    @testset "the device walks the same trajectory as the SciML loop" begin
        let tspan = (0.0, 1.0e-6), dt = 1.0e-9, saveat = 1.0e-8,
            param = prepare(zero_E, uniform_B; species = Proton),
            prob = TraceProblem([0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan, param)

            for alg in (Boris(), MultistepBoris2(n = 2), MultistepBoris4(n = 2),
                MultistepBoris6(n = 4))
                sol_loop = TestParticle.solve(prob, alg; dt, saveat)
                sol_device = TP.solve(prob, alg, CPU(); dt, trajectories = 1, saveat).u[1]

                @test length(sol_device.u) == length(sol_loop.u)
                @test sol_device.t ≈ sol_loop.t rtol = 1.0e-10
                # Both paths now call the same step, so they agree to rounding
                # rather than merely to the order of the method.
                @test max_rel_diff(sol_device.u, sol_loop.u) < 1.0e-12
            end
        end
    end

    @testset "a solver that picks its own step has no device path" begin
        let tspan = (0.0, 1.0e-6), dt = 1.0e-9,
            param = prepare(zero_E, uniform_B; species = Proton),
            prob = TraceProblem([0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan, param)

            for alg in (AdaptiveBoris(safety = 0.1), AdaptiveMultistepBoris{2}(n = 2))
                @test_throws ArgumentError TP.solve(prob, alg, CPU(); dt, trajectories = 1)
            end
        end
    end
end

end
