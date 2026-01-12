using TestParticle
using StaticArrays
using Test

@testset "Guiding Center Dispatch" begin
    # Setup
    E_field(x, t) = SA[0.0, 0.0, 0.0]
    # Time-dependent B field to verify t is correctly handled
    B_field_t(x, t) = SA[0.0, 0.0, 1.0 + t]

    param_t = prepare(E_field, B_field_t)

    # Case 1: 7-element vector (x, y, z, vx, vy, vz, t)
    # t = 1.0 -> B = 2.0
    # v = (1, 0, 0), B = (0, 0, 2)
    # B x v = (0, 2, 0)
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
