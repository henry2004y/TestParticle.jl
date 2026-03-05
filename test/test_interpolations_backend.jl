module test_interpolations_backend

using Test
using TestParticle
using StaticArrays
import TestParticle as TP

@testset "InterpolationsBackend" begin
    backend = InterpolationsBackend()

    @testset "Cartesian 3D scalar" begin
        x = range(-10, 10, length = 4)
        y = range(-10, 10, length = 6)
        z = range(-10, 10, length = 8)
        n = Float64[
            i + j + k
                for i in eachindex(x), j in eachindex(y), k in eachindex(z)
        ]

        # Periodic BC: out-of-domain point should return a finite wrapped value
        nfunc = TP.build_interpolator(backend, TP.CartesianGrid, n, x, y, z, 1, 2)
        @test isfinite(nfunc(SA[20, 0, 0]))

        # Flat BC: out-of-domain point should clamp to boundary value
        nfunc_flat = TP.build_interpolator(backend, TP.CartesianGrid, n, x, y, z, 1, 3)
        @test nfunc_flat(SA[20, 0, 0]) ≈ 12.0

        # Interior point: exact for linear interpolation
        nfunc_int = TP.build_interpolator(backend, TP.CartesianGrid, n, x, y, z, 1, 3)
        @test nfunc_int(SA[0, 0, 0]) ≈ 10.5 atol = 0.5

        # Higher orders: just verify they run and return finite values
        nfunc2 = TP.build_interpolator(backend, TP.CartesianGrid, n, x, y, z, 2, 3)
        @test isfinite(nfunc2(SA[0, 0, 0]))

        nfunc3 = TP.build_interpolator(backend, TP.CartesianGrid, n, x, y, z, 3, 3)
        @test isfinite(nfunc3(SA[0, 0, 0]))
    end

    @testset "Cartesian 3D vector" begin
        x = range(-10, 10, length = 5)
        y = range(-10, 10, length = 5)
        z = range(-10, 10, length = 5)
        B = fill(0.0, 3, length(x), length(y), length(z))
        B[3, :, :, :] .= 1.0

        Bfunc = TP.build_interpolator(backend, TP.CartesianGrid, B, x, y, z)
        @test Bfunc(SA[0.0, 0.0, 0.0]) ≈ SA[0.0, 0.0, 1.0]
    end

    @testset "RectilinearGrid" begin
        x = [0.0, 1.0, 4.0, 9.0]
        y = [0.0, 1.0, 4.0, 9.0]
        z = [0.0, 1.0, 4.0, 9.0]
        A = [i + j + k for i in x, j in y, k in z]

        Afunc = TP.build_interpolator(backend, TP.RectilinearGrid, A, x, y, z)
        @test Afunc(SA[4.0, 4.0, 4.0]) ≈ 12.0
        @test Afunc(SA[1.0, 4.0, 9.0]) ≈ 14.0
    end

    @testset "Spherical 3D" begin
        r = range(0.1, 10, length = 11)
        θ = range(0, π, length = 11)
        ϕ = range(0, 2π, length = 11)
        B = fill(0.0, 3, length(r), length(θ), length(ϕ))
        B[1, :, :, :] .= 1.0

        Bfunc = TP.build_interpolator(backend, TP.StructuredGrid, B, r, θ, ϕ)
        @test Bfunc(SA[1, 1, 1]) ≈ [0.57735, 0.57735, 0.57735] atol = 1.0e-5

        A = ones(length(r), length(θ), length(ϕ))
        Afunc = TP.build_interpolator(backend, TP.StructuredGrid, A, r, θ, ϕ)
        @test Afunc(SA[1, 1, 1]) == 1.0
    end

    @testset "prepare kwarg" begin
        x = range(-5, 5, length = 10)
        y = range(-5, 5, length = 10)
        z = range(-5, 5, length = 10)
        B = fill(0.0, 3, length(x), length(y), length(z))
        B[3, :, :, :] .= 1.0
        E = fill(0.0, 3, length(x), length(y), length(z))

        param = prepare(x, y, z, E, B; backend = InterpolationsBackend())
        @test param[4](SA[0.0, 0.0, 0.0]) ≈ SA[0.0, 0.0, 1.0]
    end

    @testset "Results match FastInterpolations" begin
        x = range(-5, 5, length = 10)
        y = range(-5, 5, length = 10)
        z = range(-5, 5, length = 10)
        n = Float64[
            sin(i) + cos(j) + i * j * k
                for i in x, j in y, k in z
        ]

        pt = SA[1.0, -1.0, 2.0]
        fi_val = TP.build_interpolator(TP.CartesianGrid, n, x, y, z)(pt)
        int_val = TP.build_interpolator(backend, TP.CartesianGrid, n, x, y, z)(pt)
        @test fi_val ≈ int_val atol = 1.0e-4
    end
end

end # module test_interpolations_backend
