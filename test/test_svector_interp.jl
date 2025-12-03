using TestParticle
using StaticArrays
using Test
using Random

@testset "SVector Interpolation" begin
    # Cartesian Grid
    nx, ny, nz = 10, 10, 10
    gridx = range(0.0, 10.0, length=nx)
    gridy = range(0.0, 10.0, length=ny)
    gridz = range(0.0, 10.0, length=nz)

    # Create random field
    Random.seed!(42)
    B_array = rand(3, nx, ny, nz)
    B_svector = reinterpret(reshape, SVector{3, Float64}, B_array)

    # Interpolator using Array{Float64, 4} (Implicitly converted internally)
    itp_implicit = TestParticle.get_interpolator(TestParticle.Cartesian(), B_array, gridx, gridy, gridz)

    # Interpolator using Array{SVector, 3} (Explicitly passed)
    itp_explicit = TestParticle.getinterp(TestParticle.Cartesian(), B_svector, gridx, gridy, gridz)

    # Compare
    pt = SA[5.5, 5.5, 5.5]
    @test itp_implicit(pt) ≈ itp_explicit(pt)
    @test itp_explicit(pt) isa SVector{3, Float64}

    # Spherical Grid
    r = range(1.0, 10.0, length=10)
    theta = range(0.0, π, length=10)
    phi = range(0.0, 2π, length=10)

    B_array_sph = rand(3, 10, 10, 10)
    B_svector_sph = reinterpret(reshape, SVector{3, Float64}, B_array_sph)

    itp_sph_implicit = TestParticle.get_interpolator(TestParticle.Spherical(), B_array_sph, r, theta, phi)
    itp_sph_explicit = TestParticle.getinterp(TestParticle.Spherical(), B_svector_sph, r, theta, phi)

    pt_cart = SA[5.0, 0.0, 0.0] # On x-axis

    @test itp_sph_implicit(pt_cart) ≈ itp_sph_explicit(pt_cart)
    @test itp_sph_explicit(pt_cart) isa SVector{3, Float64}

    # 2D case
    B_array_2d = rand(3, nx, ny)
    B_svector_2d = reinterpret(reshape, SVector{3, Float64}, B_array_2d)

    itp_2d_implicit = TestParticle.getinterp(TestParticle.Cartesian(), B_array_2d, gridx, gridy)
    itp_2d_explicit = TestParticle.getinterp(TestParticle.Cartesian(), B_svector_2d, gridx, gridy)

    pt_2d = SA[5.5, 5.5]
    @test itp_2d_implicit(pt_2d) ≈ itp_2d_explicit(pt_2d)
end
