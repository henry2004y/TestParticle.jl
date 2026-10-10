if !isdefined(Main, :test_common)
    include("test_common.jl")
end

module test_spherical_gpu

using Test
using TestParticle, KernelAbstractions, StaticArrays
import TestParticle as TP
using LinearAlgebra: norm
using ..test_common: max_rel_diff

"""
    uniform_z_field(r, θ, ϕ; B₀ = 1.0e-8) -> B

A uniform field `B₀ ẑ` expressed in spherical components, i.e.
`Br = B₀ cosθ`, `Bθ = -B₀ sinθ`, `Bϕ = 0`.
"""
function uniform_z_field(r, θ, ϕ; B₀ = 1.0e-8)
    B = zeros(3, length(r), length(θ), length(ϕ))
    for (iθ, θv) in enumerate(θ)
        s, c = sincos(θv)
        B[1, :, iθ, :] .= B₀ * c
        B[2, :, iθ, :] .= -B₀ * s
    end
    return B
end

"Scalar field that grows linearly with radius."
radial_scalar_field(r, θ, ϕ) = [1.0e-8 * rv for rv in r, _ in θ, _ in ϕ]

"""
    radial_vector_field(r, θ, ϕ; E₀ = 1.0e-8) -> E

An outward radial field `E₀ r̂` in spherical components, i.e. `Er = E₀`, `Eθ = Eϕ = 0`.
"""
function radial_vector_field(r, θ, ϕ; E₀ = 1.0e-8)
    E = zeros(3, length(r), length(θ), length(ϕ))
    E[1, :, :, :] .= E₀
    return E
end

"""
    quad_points(r, θ, ϕ)

Locations strictly inside the grid cells, at the quarter points of each direction.
"""
function quad_points(r, θ, ϕ)
    pts = SVector{3, Float64}[]
    for ir in 1:(length(r) - 1), iθ in 1:(length(θ) - 1), iϕ in 1:(length(ϕ) - 1)
        for w in (0.25, 0.5, 0.75)
            rv = r[ir] + w * (r[ir + 1] - r[ir])
            θv = θ[iθ] + w * (θ[iθ + 1] - θ[iθ])
            ϕv = ϕ[iϕ] + w * (ϕ[iϕ + 1] - ϕ[iϕ])
            push!(pts, TP.sph2cart(rv, θv, ϕv))
        end
    end
    return pts
end

"Largest deviation of the device grid from the host interpolator, relative to the field."
function max_field_rel_diff(g, fi, pts)
    return maximum(pts) do x
        a = g(x[1], x[2], x[3])
        b = fi(x)
        return norm(a - b) / max(norm(b), 1.0e-30)
    end
end

"Spherical gridded field prepared for tracing, together with its host interpolator."
function setup_field(r, θ, ϕ)
    B = uniform_z_field(r, θ, ϕ)
    A = radial_scalar_field(r, θ, ϕ)
    param = prepare(r, θ, ϕ, A, B; species = Proton, gridtype = TP.StructuredGrid)
    return param, param[4].field_function, param[3].field_function
end

@testset "GPUSphericalGrid" begin
    @testset "uniform grid agrees with the host interpolator" begin
        let r = 1.0:1.0:5.0, θ = range(0.0, π, length = 9), ϕ = range(0.0, 2π, length = 9)
            _, fiB, fiA = setup_field(r, θ, ϕ)
            gB = TP._to_gpu_spherical_grid(fiB.itp, CPU())
            gA = TP._to_gpu_spherical_grid(fiA.itp, CPU())

            @test gB isa TP.GPUSphericalGrid
            @test gB.axis_r isa TP.GPUUniformAxis
            @test gB.axis_θ isa TP.GPUUniformAxis
            @test gB.axis_ϕ isa TP.GPUUniformAxis
            @test gB.axis_r.n == Int32(length(r))
            @test gB.axis_r.x0 == first(r) && gB.axis_r.xmax == last(r)

            pts = quad_points(r, θ, ϕ)
            @test max_field_rel_diff(gB, fiB, pts) < 1.0e-12
            @test maximum(pts) do x
                return abs(gA(x[1], x[2], x[3]) - fiA(x)) / fiA(x)
            end < 1.0e-12
            # A uniform field expressed in spherical components transforms back to B₀ ẑ
            @test maximum(pts) do x
                return norm(gB(x[1], x[2], x[3]) - SA[0.0, 0.0, 1.0e-8]) / 1.0e-8
            end < 0.05
            # A field linear in r is reproduced exactly by linear interpolation
            @test maximum(pts) do x
                return abs(gA(x[1], x[2], x[3]) - 1.0e-8 * norm(x)) / (1.0e-8 * norm(x))
            end < 1.0e-12
        end
    end

    @testset "non-uniform grid agrees with the host interpolator" begin
        let r = TP.logrange(1.0, 10.0, 9), θ = range(0.0, π, length = 9),
                ϕ = range(0.0, 2π, length = 9)

            _, fiB, fiA = setup_field(r, θ, ϕ)
            gB = TP._to_gpu_spherical_grid(fiB.itp, CPU())
            gA = TP._to_gpu_spherical_grid(fiA.itp, CPU())

            @test gB.axis_r isa TP.GPUNonUniformAxis
            @test gB.axis_r.n == Int32(length(r))
            @test gB.axis_r.x0 == first(r) && gB.axis_r.xmax == last(r)

            pts = quad_points(collect(r), θ, ϕ)
            @test max_field_rel_diff(gB, fiB, pts) < 1.0e-12
            @test maximum(pts) do x
                return abs(gA(x[1], x[2], x[3]) - fiA(x)) / fiA(x)
            end < 1.0e-12
            # Bracketing a non-uniform axis must stay inside the grid
            @test all(pts) do x
                return !TP._axis_eval(gB.axis_r, norm(x))[3]
            end
        end

        # Grids handed over as plain vectors follow the same path
        let r = collect(TP.logrange(1.0, 10.0, 9)), θ = collect(range(0.0, π, length = 9)),
                ϕ = collect(range(0.0, 2π, length = 9))

            param, fiB, _ = setup_field(r, θ, ϕ)
            gB = TP._to_gpu_spherical_grid(fiB.itp, CPU())

            @test gB.axis_r isa TP.GPUNonUniformAxis
            @test max_field_rel_diff(gB, fiB, quad_points(r, θ, ϕ)) < 1.0e-12
        end

        # A non-uniform r combined with a uniform θ and ϕ keeps both axis kinds
        let r = logrange(1.0, 10.0, 16), θ = range(0.0, π, length = 16),
                ϕ = range(0.0, 2π, length = 16)

            _, fiB, _ = setup_field(r, θ, ϕ)
            gB = TP._to_gpu_spherical_grid(fiB.itp, CPU())

            @test gB.axis_r isa TP.GPUNonUniformAxis
            @test gB.axis_θ isa TP.GPUUniformAxis
            @test gB.axis_ϕ isa TP.GPUUniformAxis
            @test max_field_rel_diff(gB, fiB, quad_points(collect(r), θ, ϕ)) < 1.0e-12
        end
    end

    @testset "grid nodes are reproduced" begin
        let r = 1.0:1.0:5.0, θ = range(0.0, π, length = 9), ϕ = range(0.0, 2π, length = 9)
            B = uniform_z_field(r, θ, ϕ)
            param = prepare(
                r, θ, ϕ, zeros(size(B)), B; species = Proton, gridtype = TP.StructuredGrid
            )
            fi = param[4].field_function
            g = TP._to_gpu_spherical_grid(fi.itp, CPU())

            err = 0.0
            for (ir, rv) in enumerate(r), (iθ, θv) in enumerate(θ), (iϕ, ϕv) in enumerate(ϕ)
                x = TP.sph2cart(rv, θv, ϕv)
                Br, Bθ, Bϕ = B[1, ir, iθ, iϕ], B[2, ir, iθ, iϕ], B[3, ir, iθ, iϕ]
                expected = TP.sph2cartvec(Br, Bθ, Bϕ, θv, ϕv)
                err = max(err, norm(g(x[1], x[2], x[3]) - expected) / 1.0e-8)
            end
            # The boundary nodes sit on the domain edge, where the Cartesian to
            # spherical round trip may land an ulp outside the grid
            @test err < 1.0e-12
        end
    end

    @testset "call signatures" begin
        let r = 1.0:1.0:5.0, θ = range(0.0, π, length = 9), ϕ = range(0.0, 2π, length = 9)
            _, fiB, _ = setup_field(r, θ, ϕ)
            g = TP._to_gpu_spherical_grid(fiB.itp, CPU())
            x = SA[1.5, 2.0, 1.0]
            expected = fiB(x)

            @test g(x[1], x[2], x[3]) ≈ expected
            @test g(x[1], x[2], x[3], 0.0) ≈ expected
            @test g((x[1], x[2], x[3])) ≈ expected
            @test g((x[1], x[2], x[3]), 0.0) ≈ expected
            @test g(x) ≈ expected
            @test g(x, 0.0) ≈ expected
            @test g([x[1], x[2], x[3]]) ≈ expected

            # Adapting to a backend keeps the grid usable
            g_adapted = TP.Adapt.adapt(CPU(), g)
            @test g_adapted(x[1], x[2], x[3]) == g(x[1], x[2], x[3])
        end
    end

    @testset "outside the domain" begin
        let r = 1.0:1.0:5.0, θ = range(0.0, π, length = 9), ϕ = range(0.0, 2π, length = 9)
            _, fiB, fiA = setup_field(r, θ, ϕ)
            gB = TP._to_gpu_spherical_grid(fiB.itp, CPU())
            gA = TP._to_gpu_spherical_grid(fiA.itp, CPU())

            # Inside the radial range
            @test !any(isnan, gB(1.0, 1.0, 1.0))
            # Below r_min, above r_max, at the origin and on a NaN coordinate
            @test all(isnan, gB(0.5, 0.0, 0.0))
            @test all(isnan, gB(6.0, 0.0, 0.0))
            @test all(isnan, gB(0.0, 0.0, 0.0))
            @test all(isnan, gB(NaN, 0.0, 0.0))
            # The poles are inside the θ range
            @test !any(isnan, gB(0.0, 0.0, 2.0))
            @test !any(isnan, gB(0.0, 0.0, -2.0))
            # Scalar fields fill with a scalar NaN instead of a vector
            @test isnan(gA(0.5, 0.0, 0.0))
            @test isnan(gA(6.0, 0.0, 0.0))
            @test gA(1.0, 1.0, 1.0) ≈ 1.0e-8 * norm(SA[1.0, 1.0, 1.0])
        end
    end

    @testset "boundary conditions" begin
        # Clamping r holds the edge value instead of filling
        let r = 1.0:1.0:5.0, θ = range(0.0, π, length = 9), ϕ = range(0.0, 2π, length = 9)
            data = [SVector{3, Float64}(rv, 0.0, 0.0) for rv in r, _ in θ, _ in ϕ]
            bc_θ = TP.FillExtrap(NaN)
            ax_r = TP.GPUUniformAxis(1.0, 1.0, 1.0, Int32(length(r)), TP.ClampExtrap())
            ax_θ = TP.GPUUniformAxis(0.0, π / 8, 8 / π, Int32(length(θ)), bc_θ)
            ax_ϕ = TP.GPUUniformAxis(0.0, π / 4, 4 / π, Int32(length(ϕ)), TP.WrapExtrap())
            bc = (TP.ClampExtrap(), bc_θ, TP.WrapExtrap())
            g = TP.GPUSphericalGrid(data, ax_r, ax_θ, ax_ϕ, bc)

            θq, ϕq = π / 3, 1.0
            @test g(Tuple(TP.sph2cart(0.5, θq, ϕq))) ≈
                TP.sph2cartvec(1.0, 0.0, 0.0, θq, ϕq)
            @test g(Tuple(TP.sph2cart(8.0, θq, ϕq))) ≈
                TP.sph2cartvec(5.0, 0.0, 0.0, θq, ϕq)
            @test g(Tuple(TP.sph2cart(2.0, θq, ϕq))) ≈
                TP.sph2cartvec(2.0, 0.0, 0.0, θq, ϕq)
        end

        # Wrapping makes an axis periodic over [x0, x0 + span)
        let ax = TP.GPUUniformAxis(0.0, 1.0, 1.0, Int32(5), TP.WrapExtrap())
            i1, w1, out1 = TP._axis_eval(ax, 0.25)
            i2, w2, out2 = TP._axis_eval(ax, 4.25)
            i3, w3, out3 = TP._axis_eval(ax, -3.75)

            @test !out1 && !out2 && !out3
            @test (i1, w1) == (i2, w2) == (i3, w3)
        end

        # Filling reports the configured value in the field type
        let bc = TP.FillExtrap(NaN)
            @test all(isnan, TP._fill_value(bc, SVector{3, Float64}))
            @test isnan(TP._fill_value(bc, Float64))
        end
    end

    @testset "non-uniform axis search" begin
        let coords = [1.0, 2.0, 4.0, 8.0],
                ax = TP.GPUNonUniformAxis(coords, TP.FillExtrap(NaN))
            @test TP._axis_eval(ax, 0.5) == (Int32(1), 0.0, true)
            @test TP._axis_eval(ax, 1.0) == (Int32(1), 0.0, false)
            @test TP._axis_eval(ax, 1.5) == (Int32(1), 0.5, false)
            @test TP._axis_eval(ax, 3.0) == (Int32(2), 0.5, false)
            @test TP._axis_eval(ax, 4.0) == (Int32(3), 0.0, false)
            @test TP._axis_eval(ax, 8.0) == (Int32(3), 1.0, false)
            @test TP._axis_eval(ax, 9.0) == (Int32(1), 0.0, true)
        end
    end

    @testset "invalid input" begin
        @test_throws ArgumentError TP._to_gpu_axis(1.0, CPU(), TP.FillExtrap(NaN))
        @test_throws ArgumentError TP._to_gpu_spherical_grid((grid = 1,), CPU())
    end

    @testset "tracing on a backend" begin
        # A 1 m/s proton in a 10 nT field gyrates within a radius of about 1 m,
        # so it stays well inside the radial range of the grid.
        let r = 1.0:1.0:10.0, θ = range(0.0, π, length = 9), ϕ = range(0.0, 2π, length = 9),
                dt = 1.0e-5, tspan = (0.0, 1.0e-3)

            B = uniform_z_field(r, θ, ϕ)
            param = prepare(
                r, θ, ϕ, TP.ZeroField(), B; species = Proton, gridtype = TP.StructuredGrid
            )
            gB = TP._to_gpu_spherical_grid(param[4].field_function.itp, CPU())
            stateinit = [2.0, 2.0, 2.0, 1.0, 0.0, 0.0]

            prob_host = TraceProblem(stateinit, tspan, param)
            prob_device = TraceProblem(
                stateinit, tspan,
                (param[1], param[2], param[3], TP.Field(gB), param[5])
            )

            sol_host = TP.solve(prob_host, Boris(); dt)
            sol_device = TP.solve(prob_device, Boris(), CPU(); dt, trajectories = 1).u[1]

            @test length(sol_device.u) == length(sol_host.u)
            @test max_rel_diff(sol_device.u, sol_host.u) < 1.0e-10
        end
    end

    @testset "tracing in a spherical EM field" begin
        # Both fields live on the same spherical grid; E points radially outward,
        # B along z. A 1 m/s proton stays well inside the radial range.
        let r = 1.0:1.0:10.0, θ = range(0.0, π, length = 9), ϕ = range(0.0, 2π, length = 9),
                dt = 1.0e-4, tspan = (0.0, 0.5)

            E = radial_vector_field(r, θ, ϕ)
            B = uniform_z_field(r, θ, ϕ)
            param = prepare(r, θ, ϕ, E, B; species = Proton, gridtype = TP.StructuredGrid)
            stateinit = [2.0, 2.0, 2.0, 1.0, 0.0, 0.0]

            gE = TP._to_gpu_spherical_grid(param[3].field_function.itp, CPU())
            gB = TP._to_gpu_spherical_grid(param[4].field_function.itp, CPU())

            prob_host = TraceProblem(stateinit, tspan, param)
            prob_device = TraceProblem(
                stateinit, tspan,
                (param[1], param[2], TP.Field(gE), TP.Field(gB), param[5])
            )

            sol_host = TP.solve(prob_host, Boris(); dt)
            sol_device = TP.solve(prob_device, Boris(), CPU(); dt, trajectories = 1).u[1]

            @test length(sol_device.u) == length(sol_host.u)
            @test max_rel_diff(sol_device.u, sol_host.u) < 1.0e-10
            # The electric field does work, so the run is not a pure gyration
            @test norm(sol_device.u[end][4:6]) > norm(sol_device.u[1][4:6])
        end
    end
end

end # module test_spherical_gpu
