using TestParticle
using StaticArrays
using Test
using LinearAlgebra

@testset "GC <-> Full Conversion" begin
    # Setup simple uniform field
    B0 = 1.0 # T
    B_field(r) = SA[0, 0, B0]
    E_field(r) = SA[0, 0, 0]
    param = prepare(E_field, B_field; species = Proton)

    @testset "Reversibility GC -> Full -> GC" begin
        # Define GC state
        # R = [1.0, 0.0, 0.0]
        # vpar = 1.0e5
        # mu = 1.0e-5 # some value

        R = SA[1.0, 0.0, 0.0]
        vpar = 1.0e5
        gc_state = [R..., vpar]
        mu = 1.0e-15

        # Random phase
        phase = π / 3

        # Convert to full
        xu = gc_to_full(gc_state, param, mu, phase)

        # Check dimensions
        @test length(xu) == 6

        # Convert back to GC
        gc_new, mu_new = full_to_gc(xu, param)

        # Check consistency
        @test gc_new[1:3] ≈ R rtol = 1.0e-5
        @test gc_new[4] ≈ vpar rtol = 1.0e-5
        @test mu_new ≈ mu rtol = 1.0e-5
    end

    @testset "With E-field (ExB drift)" begin
        # B in Z, E in Y -> ExB in X
        B0 = 1.0
        E0 = 1.0e3
        B_field2(r) = SA[0, 0, B0]
        E_field2(r) = SA[0, E0, 0] # E in y
        param2 = prepare(E_field2, B_field2; species = Proton)

        # ExB drift: v_E = (E x B) / B^2
        # (0, E0, 0) x (0, 0, B0) = (E0*B0, 0, 0) -> direction x
        # magnitude E0*B0 / B0^2 = E0/B0 = 1000.

        R = SA[0.0, 0.0, 0.0]
        vpar = 0.0
        gc_state = [R..., vpar]
        mu = 0.0 # No gyromotion, only drift

        xu = gc_to_full(gc_state, param2, mu)

        # v should be just drift velocity + vpar(0)
        v = xu[4:6]
        @test v[1] ≈ 1000.0 atol = 1.0e-5
        @test v[2] ≈ 0.0 atol = 1.0e-5
        @test v[3] ≈ 0.0 atol = 1.0e-5

        # Reverse check
        gc_new, mu_new = full_to_gc(xu, param2)
        @test gc_new[1:3] ≈ R atol = 1.0e-5 # Position should match GC
        @test gc_new[4] ≈ vpar atol = 1.0e-5
        @test mu_new ≈ 0.0 atol = 1.0e-10
    end

    @testset "Phase check" begin
        # Check that phase argument changes velocity direction
        B_field3(r) = SA[0, 0, 1.0]
        param3 = prepare(E_field, B_field3; species = Proton)

        R = SA[0, 0, 0]
        vpar = 0.0
        gc_state = [R..., vpar]
        mu = 1.0e-15

        # Phase 0
        xu0 = gc_to_full(gc_state, param3, mu, 0.0)
        v0 = xu0[4:6]

        # Phase pi/2
        xu90 = gc_to_full(gc_state, param3, mu, π / 2)
        v90 = xu90[4:6]

        # For B in Z, perp plane is XY.
        # e1, e2 from get_perp_vector([0,0,1]):
        # v = [0,1,0], e1 = v x b = [1,0,0], e2 = b x e1 = [0,1,0].
        # So phase 0 -> x direction. Phase 90 -> y direction.


        # Check orthogonality
        # Normalize to avoid scaling issues with large w
        @test dot(normalize(v0), normalize(v90)) ≈ 0.0 atol = 1.0e-5

        # Verify e1 aligned with x (since b=z, v_arbitrary=y -> e1=x)
        # v0 should be in x direction
        @test isapprox(normalize(v0), SA[1.0, 0.0, 0.0], atol = 1.0e-5) ||
            isapprox(normalize(v0), SA[-1.0, 0.0, 0.0], atol = 1.0e-5)

        # v90 should be in y direction
        @test isapprox(normalize(v90), SA[0.0, 1.0, 0.0], atol = 1.0e-5) ||
            isapprox(normalize(v90), SA[0.0, -1.0, 0.0], atol = 1.0e-5)
    end
end
