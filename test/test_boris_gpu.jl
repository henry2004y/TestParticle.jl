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
                sol_loop = solve(prob, alg; dt, saveat)
                sol_device = solve(prob, alg, CPU(); dt, trajectories = 1, saveat).u[1]

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
                @test_throws ArgumentError solve(prob, alg, CPU(); dt, trajectories = 1)
            end
        end
    end

    @testset "Float32 Boris on a backend" begin
        let tspan = (0.0f0, 1.0f-6), dt = 1.0f-9, saveat = 1.0f-8,
                param = prepare(zero_E, uniform_B; species = Proton),
                prob = TraceProblem(Float32[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0], tspan, param)

            sol_loop = solve(prob, Boris(); dt, saveat)
            sol_device = solve(prob, Boris(), CPU(); dt, trajectories = 1, saveat).u[1]

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

            sol_loop = solve(prob, Boris(); dt = 1.0e-9)
            sol_device = solve(prob, Boris(), CPU(); dt = 1.0e-9, trajectories = 1).u[1]

            @test length(sol_device.u) == length(sol_loop.u)
            @test max_rel_diff(sol_device.u, sol_loop.u) < 1.0e-12

            if isdefined(Main, :CUDA) && Main.CUDA.functional()
                CUDA = Main.CUDA
                backend = CUDA.CUDABackend()
                sol_gpu = solve(
                    prob, Boris(), backend; dt = 1.0e-9, trajectories = 1
                ).u[1]
                @test length(sol_gpu.u) == length(sol_loop.u)
                @test max_rel_diff(sol_gpu.u, sol_loop.u) < 1.0e-12

                # Pinned host staging buffer
                ens_pin = EnsembleKernel(backend; pin = true)
                @test ens_pin.pin == true
                sol_gpu_pin = solve(
                    prob, Boris(), ens_pin;
                    dt = 1.0e-9, trajectories = 2, raw_output = true
                )
                @test size(sol_gpu_pin.u, 1) == 2
                @test length(sol_gpu_pin.t) == length(sol_gpu.t)
            end
        end

        # Spherical Grid (Uniform and Non-Uniform)
        let r_u = range(1.0, 10.0, length = 11),
                θ = range(0.0, π, length = 11),
                ϕ = range(0.0, 2π, length = 11),
                B_sph = zeros(3, length(r_u), length(θ), length(ϕ))

            for (iθ, θ_val) in enumerate(θ)
                sinθ, cosθ = sincos(θ_val)
                B_sph[1, :, iθ, :] .= 1.0e-8 * cosθ
                B_sph[2, :, iθ, :] .= -1.0e-8 * sinθ
            end

            # Uniform spherical grid
            param_sph_u = prepare(
                r_u, θ, ϕ, zeros(size(B_sph)), B_sph;
                species = Proton, gridtype = TP.StructuredGrid
            )
            prob_sph_u = TraceProblem(
                [2.0, 2.0, 2.0, 1.0, 0.0, 0.0], (0.0, 1.0e-3), param_sph_u
            )
            sol_loop_u = solve(prob_sph_u, Boris(); dt = 1.0e-5)
            sol_cpu_u = solve(
                prob_sph_u, Boris(), CPU(); dt = 1.0e-5, trajectories = 1
            ).u[1]

            @test length(sol_cpu_u.u) == length(sol_loop_u.u)
            @test max_rel_diff(sol_cpu_u.u, sol_loop_u.u) < 1.0e-12

            # Non-uniform spherical grid (LogRange in r)
            r_nu = TP.logrange(1.0, 10.0, 11)
            B_sph_nu = zeros(3, length(r_nu), length(θ), length(ϕ))
            for (iθ, θ_val) in enumerate(θ)
                sinθ, cosθ = sincos(θ_val)
                B_sph_nu[1, :, iθ, :] .= 1.0e-8 * cosθ
                B_sph_nu[2, :, iθ, :] .= -1.0e-8 * sinθ
            end
            param_sph_nu = prepare(
                r_nu, θ, ϕ, zeros(size(B_sph_nu)), B_sph_nu;
                species = Proton, gridtype = TP.StructuredGrid
            )
            prob_sph_nu = TraceProblem(
                [2.0, 2.0, 2.0, 1.0, 0.0, 0.0], (0.0, 1.0e-3), param_sph_nu
            )
            sol_loop_nu = solve(prob_sph_nu, Boris(); dt = 1.0e-5)
            sol_cpu_nu = solve(
                prob_sph_nu, Boris(), CPU(); dt = 1.0e-5, trajectories = 1
            ).u[1]

            @test length(sol_cpu_nu.u) == length(sol_loop_nu.u)
            @test max_rel_diff(sol_cpu_nu.u, sol_loop_nu.u) < 1.0e-12

            if isdefined(Main, :CUDA) && Main.CUDA.functional()
                CUDA = Main.CUDA
                backend = CUDA.CUDABackend()

                sol_gpu_u = solve(
                    prob_sph_u, Boris(), backend; dt = 1.0e-5, trajectories = 2
                ).u[1]
                @test length(sol_gpu_u.u) == length(sol_loop_u.u)
                @test max_rel_diff(sol_gpu_u.u, sol_loop_u.u) < 1.0e-5

                sol_gpu_nu = solve(
                    prob_sph_nu, Boris(), backend; dt = 1.0e-5, trajectories = 2
                ).u[1]
                @test length(sol_gpu_nu.u) == length(sol_loop_nu.u)
                @test max_rel_diff(sol_gpu_nu.u, sol_loop_nu.u) < 1.0e-5
            end
        end
    end

    @testset "GPUGrid interpolators on CPU" begin
        x = 0.0:1.0:4.0
        y = 0.0:1.0:4.0
        z = 0.0:1.0:4.0

        # 3D Grid
        B_grid3d = zeros(3, length(x), length(y), length(z))
        for i in 1:length(x), j in 1:length(y), k in 1:length(z)
            B_grid3d[:, i, j, k] .= [x[i], y[j], z[k]]
        end
        param3d = prepare(x, y, z, zeros(size(B_grid3d)), B_grid3d; species = Proton)
        fi3d = param3d[4].field_function
        g3d = TP._to_gpu_grid(fi3d.itp, CPU())

        @test g3d(1.5, 2.5, 3.5) ≈ SA[1.5, 2.5, 3.5]
        @test g3d((1.5, 2.5, 3.5)) ≈ SA[1.5, 2.5, 3.5]
        @test g3d(SA[1.5, 2.5, 3.5]) ≈ SA[1.5, 2.5, 3.5]
        @test g3d(SA[1.5, 2.5, 3.5], 0.0) ≈ SA[1.5, 2.5, 3.5]

        g3d_adapt = TP.Adapt.adapt(CPU(), g3d)
        @test g3d_adapt(1.0, 1.0, 1.0) == g3d(1.0, 1.0, 1.0)

        # 2D Grid
        B_grid2d = zeros(3, length(x), length(y))
        for i in 1:length(x), j in 1:length(y)
            B_grid2d[:, i, j] .= [x[i], y[j], 0.0]
        end
        param2d = prepare(x, y, zeros(size(B_grid2d)), B_grid2d; species = Proton)
        fi2d = param2d[4].field_function
        g2d = TP._to_gpu_grid(fi2d.itp, CPU())

        @test g2d(1.5, 2.5) ≈ SA[1.5, 2.5, 0.0]
        @test g2d((1.5, 2.5)) ≈ SA[1.5, 2.5, 0.0]
        @test g2d(SA[1.5, 2.5, 0.0]) ≈ SA[1.5, 2.5, 0.0]
        @test g2d(SA[1.5, 2.5, 0.0], 0.0) ≈ SA[1.5, 2.5, 0.0]

        g2d_adapt = TP.Adapt.adapt(CPU(), g2d)
        @test g2d_adapt(1.0, 1.0) == g2d(1.0, 1.0)

        # 1D Grid
        B_grid1d = zeros(3, length(x))
        for i in 1:length(x)
            B_grid1d[:, i] .= [x[i], 0.0, 0.0]
        end
        param1d = prepare(x, zeros(size(B_grid1d)), B_grid1d; species = Proton)
        fi1d = param1d[4].field_function
        g1d = TP._to_gpu_grid(fi1d.itp, CPU(); dir = 1)

        @test g1d(2.5) ≈ SA[2.5, 0.0, 0.0]
        @test g1d((2.5,)) ≈ SA[2.5, 0.0, 0.0]
        @test g1d(SA[2.5, 0.0, 0.0]) ≈ SA[2.5, 0.0, 0.0]
        @test g1d(SA[2.5, 0.0, 0.0], 0.0) ≈ SA[2.5, 0.0, 0.0]

        g1d_adapt = TP.Adapt.adapt(CPU(), g1d)
        @test g1d_adapt(2.0) == g1d(2.0)

        # Spherical Grid Interpolators
        let r = 1.0:1.0:5.0, θ = range(0.0, π, length = 5), ϕ = range(0.0, 2π, length = 5)
            B_sph = zeros(3, length(r), length(θ), length(ϕ))
            B_sph[1, :, :, :] .= 1.0
            param_sph = prepare(
                r, θ, ϕ, zeros(size(B_sph)), B_sph;
                species = Proton, gridtype = TP.StructuredGrid
            )
            fi_sph = param_sph[4].field_function
            g_sph = TP._to_gpu_spherical_grid(fi_sph.itp, CPU())

            @test g_sph(1.0, 1.0, 1.0) ≈ fi_sph(SA[1.0, 1.0, 1.0])
            @test g_sph((1.0, 1.0, 1.0)) ≈ fi_sph(SA[1.0, 1.0, 1.0])
            @test g_sph(SA[1.0, 1.0, 1.0]) ≈ fi_sph(SA[1.0, 1.0, 1.0])
            @test g_sph(SA[1.0, 1.0, 1.0], 0.0) ≈ fi_sph(SA[1.0, 1.0, 1.0])
            @test isnan(g_sph(0.0, 0.0, 0.0)[1])
            @test isnan(g_sph(NaN, 0.0, 0.0)[1])

            g_sph_adapt = TP.Adapt.adapt(CPU(), g_sph)
            @test g_sph_adapt(2.0, 2.0, 2.0) == g_sph(2.0, 2.0, 2.0)

            # Non-uniform axis search bracket test
            coords = [1.0, 2.0, 4.0, 8.0]
            ax_nu = TP.GPUNonUniformAxis(coords, TP.FillExtrap(NaN))
            idx, w, out = TP._axis_eval(ax_nu, 3.0)
            @test !out
            @test idx == 2
            @test w ≈ 0.5
        end

        # Boundary conditions
        let data = collect(1.0:5.0)
            bc_fill = TP.FillExtrap(0.0)
            grid_fill = TP.GPUGrid1D(data, 1.0, 1.0, 1.0, Int32(5), Int32(1), bc_fill)
            @test grid_fill(0.0) == 0.0
            @test grid_fill(6.0) == 0.0
            @test grid_fill(2.5) == 2.5

            bc_clamp = TP.ClampExtrap()
            grid_clamp = TP.GPUGrid1D(data, 1.0, 1.0, 1.0, Int32(5), Int32(1), bc_clamp)
            @test grid_clamp(0.0) == 1.0
            @test grid_clamp(6.0) == 5.0

            bc_wrap = TP.WrapExtrap()
            grid_wrap = TP.GPUGrid1D(data, 1.0, 1.0, 1.0, Int32(5), Int32(1), bc_wrap)
            @test grid_wrap(5.5) ≈ grid_wrap(1.5)
        end

        # Tracing with Field(g3d) on CPU
        prob3d = TraceProblem(
            [1.0, 1.0, 1.0, 1.0e5, 0.0, 0.0], (0.0, 1.0e-6),
            (param3d[1], param3d[2], param3d[3], TP.Field(g3d))
        )
        sol3d = solve(prob3d, Boris(), CPU(); dt = 1.0e-9)
        @test length(sol3d.u[1].u) == 1001
    end
end

end
