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

            for alg in (
                    Boris(), MultistepBoris2(n = 2), MultistepBoris4(n = 2),
                    MultistepBoris6(n = 4),
                )
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

    @testset "Float32 Boris on a backend" begin
        let tspan = (0.0f0, 1.0f-6), dt = 1.0f-9, saveat = 1.0f-8,
                param = prepare(zero_E, uniform_B; species = Proton),
                prob = TraceProblem(Float32[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan, param)

            sol_loop = TestParticle.solve(prob, Boris(); dt, saveat)
            sol_device = TP.solve(prob, Boris(), CPU(); dt, trajectories = 1, saveat).u[1]

            @test eltype(sol_device.u[end]) === Float32
            @test eltype(sol_device.t) === Float32
            @test max_rel_diff(sol_device.u, sol_loop.u) < 1.0e-5
        end
    end

    @testset "GPU grid interpolation" begin
        let x = 0.0:1.0:10.0, y = 0.0:1.0:10.0, z = 0.0:1.0:10.0,
                B_grid = zeros(3, length(x), length(y), length(z))

            B_grid[3, :, :, :] .= 1.0e-8
            param = prepare(x, y, z, zeros(size(B_grid)), B_grid; species = Proton)
            prob = TraceProblem([5.0, 5.0, 5.0, 1.0e5, 0.0, 0.0], (0.0, 1.0e-6), param)

            sol_loop = TestParticle.solve(prob, Boris(); dt = 1.0e-9)
            sol_device = TP.solve(prob, Boris(), CPU(); dt = 1.0e-9, trajectories = 1).u[1]

            @test length(sol_device.u) == length(sol_loop.u)
            @test max_rel_diff(sol_device.u, sol_loop.u) < 1.0e-12

            if isdefined(Main, :CUDA) && Main.CUDA.functional()
                CUDA = Main.CUDA
                backend = CUDA.CUDABackend()
                sol_gpu = TP.solve(
                    prob, Boris(), backend; dt = 1.0e-9, trajectories = 1
                ).u[1]
                @test length(sol_gpu.u) == length(sol_loop.u)
                @test max_rel_diff(sol_gpu.u, sol_loop.u) < 1.0e-12
            end
        end
    end

    @testset "Morton 3D particle sorting" begin
        let N = 100,
                xv = rand(N, 6)

            perm = morton_sort_particles(xv)
            @test length(perm) == N
            @test sort(perm) == 1:N

            u0s = [xv[i, :] for i in 1:N]
            param = prepare(zero_E, uniform_B; species = Proton)
            prob_func = (prob, i, repeat = false) -> remake(prob; u0 = u0s[i.sim_id])
            base_prob = TraceProblem(u0s[1], (0.0, 1.0e-6), param; prob_func)

            sol_unsorted = TP.solve(
                base_prob, Boris(), CPU();
                trajectories = N, dt = 1.0e-9, sort_particles = false,
            )
            sol_sorted = TP.solve(
                base_prob, Boris(), CPU();
                trajectories = N, dt = 1.0e-9, sort_particles = true,
            )

            for i in 1:N
                @test sol_unsorted.u[i].u[end] ≈ sol_sorted.u[i].u[end]
            end
        end
    end
end

end
