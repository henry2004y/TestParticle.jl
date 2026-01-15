using TestParticle
using Test
using StaticArrays
using LinearAlgebra

@testset "Loop Field" begin
    # Parameters
    R = [0.0, 0.0, 0.0]
    a = 1.0
    I = 1.0
    n = [0.0, 0.0, 1.0]

    # Test field at the center
    r_center = [0.0, 0.0, 0.0]
    B_center = TestParticle.getB_loop(r_center, R, a, I, n)
    # B = mu0 * I / (2 * a) in z direction
    B_expected = TestParticle.μ₀ * I / (2 * a) * n
    @test B_center ≈ B_expected atol = 1.0e-15

    # Test field on the axis
    z = 2.0
    r_axis = [0.0, 0.0, z]
    B_axis = TestParticle.getB_loop(r_axis, R, a, I, n)
    # B = mu0 * I * a^2 / (2 * (a^2 + z^2)^(3/2))
    B_mag = TestParticle.μ₀ * I * a^2 / (2 * (a^2 + z^2)^1.5)
    @test B_axis ≈ B_mag * n atol = 1.0e-15

    # Test arbitrary orientation
    n_rot = [1.0, 0.0, 0.0]
    B_center_rot = TestParticle.getB_loop(r_center, R, a, I, n_rot)
    @test B_center_rot ≈ TestParticle.μ₀ * I / (2 * a) * n_rot atol = 1.0e-15

    # Test off-axis field consistency (check B_z at z=0 plane)
    # For z=0, B_z = mu0*I/(2*pi*a) * K(k^2) ? No, let's use the formula.
    # At z=0, B_rho = 0.
    # We can check continuity or just that it runs without error and gives reasonable numbers.

    r_off = [0.5, 0.0, 0.0]
    B_off = TestParticle.getB_loop(r_off, R, a, I, n)
    @test B_off[1] ≈ 0.0 atol = 1.0e-15 # No radial component if z=0? Wait.
    # If z=0, B_rho component:
    # factor * (z / rho) * ... -> 0 because z=0. Correct.
    # So B should be purely in z direction.
    @test abs(B_off[3]) > 0.0

    # Test rotated loop with off-axis point
    # Rotate loop by 90 degrees around Y axis, normal becomes X.
    n_x = [1.0, 0.0, 0.0]
    # Point that was at [0, 0, 2] (z-axis) should now be at [2, 0, 0] (x-axis) relative to loop
    r_x = [2.0, 0.0, 0.0]
    B_x = TestParticle.getB_loop(r_x, R, a, I, n_x)

    # Magnitude should match the previous axis test
    @test norm(B_x) ≈ B_mag atol = 1.0e-15
    @test B_x ≈ B_mag * n_x atol = 1.0e-15

    # Test translated loop
    R_new = [10.0, 10.0, 10.0]
    r_new = R_new + [0.0, 0.0, 2.0] # 2.0 along z from center
    B_trans = TestParticle.getB_loop(r_new, R_new, a, I, n)
    @test B_trans ≈ B_mag * n atol = 1.0e-15
end
