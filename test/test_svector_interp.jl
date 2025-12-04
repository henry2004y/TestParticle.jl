import TestParticle as TP
using StaticArrays
using Test
using Random

@testset "SVector Interpolation" begin
   # Cartesian Grid
   nx, ny, nz = 10, 10, 10
   gridx = range(0.0, 10.0, length = nx)
   gridy = range(0.0, 10.0, length = ny)
   gridz = range(0.0, 10.0, length = nz)

   # Create random field
   Random.seed!(42)
   B_array = rand(3, nx, ny, nz)
   B_svector = reinterpret(reshape, SVector{3, Float64}, B_array)

   # Interpolator using Array{Float64, 4} (Implicitly converted internally)
   itp_implicit = TP.get_interpolator(TP.Cartesian(), B_array, gridx, gridy, gridz)

   # Interpolator using Array{SVector, 3} (Explicitly passed)
   itp_explicit = TP.getinterp(TP.Cartesian(), B_svector, gridx, gridy, gridz)

   # Compare
   pt = SA[5.5, 5.5, 5.5]
   @test itp_implicit(pt) ≈ itp_explicit(pt)
   @test itp_explicit(pt) isa SVector{3, Float64}

   # Spherical Grid
   r = range(1.0, 10.0, length = 10)
   theta = range(0.0, π, length = 10)
   phi = range(0.0, 2π, length = 10)

   B_array_sph = rand(3, 10, 10, 10)
   B_svector_sph = reinterpret(reshape, SVector{3, Float64}, B_array_sph)

   itp_sph_implicit = TP.get_interpolator(TP.Spherical(), B_array_sph, r, theta, phi)
   itp_sph_explicit = TP.getinterp(TP.Spherical(), B_svector_sph, r, theta, phi)

   pt_cart = SA[5.0, 0.0, 0.0] # On x-axis

   @test itp_sph_implicit(pt_cart) ≈ itp_sph_explicit(pt_cart)
   @test itp_sph_explicit(pt_cart) isa SVector{3, Float64}

   # 2D case
   B_array_2d = rand(3, nx, ny)
   B_svector_2d = reinterpret(reshape, SVector{3, Float64}, B_array_2d)

   itp_2d_implicit = TP.getinterp(TP.Cartesian(), B_array_2d, gridx, gridy)
   itp_2d_explicit = TP.getinterp(TP.Cartesian(), B_svector_2d, gridx, gridy)

   pt_2d = SA[5.5, 5.5]
   @test itp_2d_implicit(pt_2d) ≈ itp_2d_explicit(pt_2d)
end

@testset "CartesianNonUniform" begin
   x = range(0.0, 10.0, length=11)
   y = range(0.0, 10.0, length=11)
   z = range(0.0, 10.0, length=11)

   # Vector field
   B = fill(0.0, 3, length(x), length(y), length(z))
   B[1, :, :, :] .= 1.0 # Bx = 1.0

   Bfunc = TP.getinterp(TP.CartesianNonUniform(), B, x, y, z)
   @test Bfunc(SA[5.0, 5.0, 5.0]) ≈ [1.0, 0.0, 0.0]

   # Scalar field
   A = fill(2.0, length(x), length(y), length(z))
   Afunc = TP.getinterp_scalar(TP.CartesianNonUniform(), A, x, y, z)
   @test Afunc(SA[5.0, 5.0, 5.0]) ≈ 2.0

   # Check interpolation values
   # Linear gradient
   A_grad = [i + j + k for i in x, j in y, k in z]
   Afunc_grad = TP.getinterp_scalar(TP.CartesianNonUniform(), A_grad, x, y, z)

   # Center point: 5.0, 5.0, 5.0 -> should match exactly for linear interpolation
   # value = 5.0 + 5.0 + 5.0 = 15.0
   @test Afunc_grad(SA[5.0, 5.0, 5.0]) ≈ 15.0

   # Off-grid point: 5.5, 5.5, 5.5
   # Should be 5.5 + 5.5 + 5.5 = 16.5
   @test Afunc_grad(SA[5.5, 5.5, 5.5]) ≈ 16.5

   # SVector array support
   B_sv = [SA[1.0, 0.0, 0.0] for i in x, j in y, k in z]
   Bfunc_sv = TP.getinterp(TP.CartesianNonUniform(), B_sv, x, y, z)
   @test Bfunc_sv(SA[5.0, 5.0, 5.0]) ≈ SA[1.0, 0.0, 0.0]

   # Non-uniform grid check
   x_nu = [0.0, 1.0, 4.0, 9.0]
   y_nu = [0.0, 1.0, 4.0, 9.0]
   z_nu = [0.0, 1.0, 4.0, 9.0]
   A_nu = [sqrt(i) + sqrt(j) + sqrt(k) for i in x_nu, j in y_nu, k in z_nu]

   # Check if we can interpolate on this
   Afunc_nu = TP.getinterp_scalar(TP.CartesianNonUniform(), A_nu, x_nu, y_nu, z_nu)

   # Point (4.0, 4.0, 4.0) -> sqrt(4)+sqrt(4)+sqrt(4) = 2+2+2=6
   @test Afunc_nu(SA[4.0, 4.0, 4.0]) ≈ 6.0

   # Point (1.0, 4.0, 9.0) -> 1+2+3=6
   @test Afunc_nu(SA[1.0, 4.0, 9.0]) ≈ 6.0
end
