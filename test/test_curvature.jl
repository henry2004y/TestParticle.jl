using Test
using StaticArrays
using TestParticle
using LinearAlgebra

function B_circular(x, t)
    r = hypot(x[1], x[3])
    r == 0 && return SA[0.0, 0.0, 0.0]
    # Magnitude = r
    # Vector = (-z, 0, x). Norm = r.
    # B = (-z, 0, x)
    return SA[-x[3], 0.0, x[1]]
end

@testset "Curvature Radius" begin
    @testset "Zero Field" begin
        B_zero(x, t) = SA[0.0, 0.0, 0.0]
        x = SA[1.0, 0.0, 0.0]
        @test get_curvature_radius(x, 0.0, B_zero) == Inf
    end

    @testset "Uniform Field" begin
        B_uniform(x, t) = SA[0.0, 0.0, 1.0]
        x = SA[1.0, 0.0, 0.0]
        @test get_curvature_radius(x, 0.0, B_uniform) == Inf
    end

    @testset "Circular Field" begin
        # Field lines are circles in x-z plane: B = (-z, 0, x)
        # Radius of curvature at (R, 0, 0) is R.
        R = 5.0
        @test get_curvature_radius(SA[R, 0.0, 0.0], 0.0, B_circular) ≈ R
        @test get_curvature_radius(SA[0.0, 0.0, R], 0.0, B_circular) ≈ R

        theta = π / 4
        @test get_curvature_radius(
            SA[R * cos(theta), 0.0, R * sin(theta)], 0.0, B_circular
        ) ≈ R
    end

    @testset "Adiabaticity" begin
        # Adiabaticity ϵ = ρ / Rc.
        # With R=5, q=1, m=2, μ=0.5:
        # Rc = R = 5.0, |B| = R = 5.0
        # ρ = sqrt(2*μ*m/|B|) / q ≈ 0.632
        # ϵ ≈ 0.126

        R = 5.0
        r = SA[R, 0.0, 0.0]

        q = 1.0
        m = 2.0
        μ = 0.5

        val = get_adiabaticity(r, B_circular, q, m, μ, 0.0)
        ρ_expected = sqrt(2 * μ * m / R) / q
        Rc_expected = R
        @test val ≈ ρ_expected / Rc_expected

        # Test dispatch
        ts = TestParticle.Species(m, q)
        @test get_adiabaticity(r, B_circular, μ, 0.0; species = ts) ≈ val

        # Test Default Proton
        # Proton: q=qᵢ, m=mᵢ
        # val_p = get_adiabaticity(r, B_circular, 1e-19, 0.0) # μ ~ 1e-19
        # just check it runs
        @test isfinite(get_adiabaticity(r, B_circular, 1.0e-19, 0.0))

        # Test Negative Charge
        # Should result in positive adiabaticity due to abs(q)
        q_neg = -1.0
        val_neg = get_adiabaticity(r, B_circular, q_neg, m, μ, 0.0)
        @test val_neg > 0
        @test val_neg ≈ val

        # Test Zero B Field
        # Should return Inf
        B_zero_test(x, t) = SA[0.0, 0.0, 0.0]
        @test get_adiabaticity(r, B_zero_test, q, m, μ, 0.0) == Inf
    end
end
