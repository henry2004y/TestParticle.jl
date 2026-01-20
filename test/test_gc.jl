using TestParticle
using StaticArrays
using Test
using LinearAlgebra

@testset "GC" begin

    @testset "Guiding Center Dispatch" begin
        # Setup
        E_field(x, t) = SA[0.0, 0.0, 0.0]
        # Time-dependent B field to verify t is correctly handled
        B_field_t(x, t) = SA[0.0, 0.0, 1.0 + t]

        param_t = prepare(E_field, B_field_t)

        # Case 1: 7-element vector (x, y, z, vx, vy, vz, t)
        # t = 1.0 -> B = 2.0
        # v = (1, 0, 0), B = (0, 0, 2), B x v = (0, 2, 0)
        # rho = B x v / (q2m * B^2)
        # q2m approx 1e8 for proton.
        # This is just to check that t is read correctly.

        xu_7 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
        gc_7 = get_gc(xu_7, param_t)

        # Case 2: 6-element vector (x, y, z, vx, vy, vz)
        # t should default to 0.0 -> B = 1.0
        # If bug exists, t would be taken as vz = 2.0 -> B = 3.0
        xu_6 = [0.0, 0.0, 0.0, 1.0, 0.0, 2.0]
        gc_6 = get_gc(xu_6, param_t)

        # For gc_7 (t=1.0): B=2.0
        # For gc_6 (t=0.0): B=1.0. If bug, t=2.0 -> B=3.0.

        # We can check the Larmor radius or the shift.
        # rho ~ 1/B.
        # shift_7 = |gc - x| ~ 1/2
        # shift_6 = |gc - x| ~ 1/1

        # Let's check B field value indirectly by checking the shift magnitude (which is Larmor radius).
        # We don't have direct access to B inside get_gc return.

        # shift = gc - x
        # x is (0,0,0)
        shift_7 = gc_7[1:3]
        shift_6 = gc_6[1:3]

        # Expected ratio of shifts: rho_6 / rho_7 = (1/B_6) / (1/B_7) = (1/1) / (1/2) = 2.
        # If bug (t=2.0 for gc_6): rho_6_bug ~ 1/3.
        # Ratio would be (1/3) / (1/2) = 2/3.

        norm_shift_7 = sqrt(sum(shift_7 .^ 2))
        norm_shift_6 = sqrt(sum(shift_6 .^ 2))

        @test isapprox(norm_shift_6 / norm_shift_7, 2.0, rtol = 0.01)
    end

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
            @test isapprox(normalize(v0), SA[1.0, 0.0, 0.0], atol = 1.0e-5)

            # v90 should be in y direction
            @test isapprox(normalize(v90), SA[0.0, 1.0, 0.0], atol = 1.0e-5)
        end
    end
end
