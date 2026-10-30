module test_hybrid

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using TestParticle
import TestParticle as TP
using StaticArrays
using LinearAlgebra
using BorisPushers
using Test
using ..test_common: bottle_B, bottle_B0, sheared_B

const m = TP.mᵢ
const q = TP.qᵢ
const q2m = q / m

const E_field = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

"Uniform 0.01 T field along +z, flat enough to stay adiabatic."
B_flat(x, t) = SA[0.0, 0.0, 0.01]

const B_field_flat = TP.Field(B_flat)
const B_field_bottle = TP.Field(bottle_B)
const B_field_sheared = TP.Field(sheared_B)

"Gyroperiod in a field of magnitude `B0`."
gyroperiod(B0) = 2π / abs(q2m * B0)

"Kinetic energy [J] of the velocity part of a full-orbit state."
kinetic_energy(v) = 0.5 * m * sum(abs2, v)

@testset "Adaptive Hybrid Solver" begin
    # 1. Pure Adiabatic Case (Uniform Field)
    @testset "uniform field" begin
        u0 = vcat(SA[0.0, 0.0, 0.0], SA[1.0e4, 0.0, 1.0e3])
        tspan = (0.0, 1.0e-4)
        p = (q2m, m, E_field, B_field_flat, TP.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)

        sols = solve(TraceProblem(u0, tspan, p), alg)
        @test sols.u[1].retcode == TP.ReturnCode.Success
    end

    # 2. Non-adiabatic Case (Highly Curved Field)
    @testset "sheared field" begin
        u0 = vcat(SA[0.0, 0.0, 0.0], SA[1.0e4, 0.0, 1.0e3])
        tspan = (0.0, 1.0e-5)
        p = (q2m, m, E_field, B_field_sheared, TP.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)

        sols = solve(TraceProblem(u0, tspan, p), alg)
        @test sols.u[1].retcode == TP.ReturnCode.Success
    end

    # 3. Demo Hybrid Case (Explicit Switching)
    # Replicates the setup from docs/examples/features/demo_hybrid.jl; the
    # stronger curvature (α = 1.0e-2) forces GC <-> FO switching.
    @testset "magnetic bottle" begin
        v0 = SA[5.0e4, 0.0, 1.0e5]
        u0 = vcat(SA[0.0, 0.0, 0.0], v0)
        T_gyro = gyroperiod(bottle_B0)
        tspan = (0.0, 30 * T_gyro)
        p = (q2m, m, E_field, B_field_bottle, TP.ZeroField())

        # Consistent with demo_hybrid.jl
        alg = AdaptiveHybrid(;
            threshold = 0.1,
            dtmax = T_gyro,
            dtmin = 1.0e-4 * T_gyro,
            check_interval = 100,
        )

        sols = solve(TraceProblem(u0, tspan, p), alg)

        @test sols.u[1].retcode == TP.ReturnCode.Success
        # Check that we actually have a reasonable number of steps (hybrid should adapt)
        @test length(sols.u[1].t) > 100

        # Energy conservation check (approximate for hybrid)
        KE_final = kinetic_energy(sols.u[1].u[end][SA[4, 5, 6]])
        @test isapprox(KE_final, kinetic_energy(v0); rtol = 0.1)
    end

    # 4. Round-trip: GC -> full-orbit reconstruction is exact in a uniform field.
    # (The guiding-center map is the exact inverse of itself when ∇B = 0.)
    @testset "GC round trip" begin
        for x0 in (SA[0.0, 0.0, 0.0], SA[1.0, 0.5, 1.0], SA[3.0, 0.0, 2.0])
            for v0 in (SA[5.0e4, 0.0, 1.0e5], SA[3.0e4, 4.0e4, 8.0e4])
                xv = vcat(x0, v0)
                X, vpar, _, _, μ, phase = TP._get_gc_parameters(
                    xv, E_field, B_field_flat, q, m, 0.0
                )
                state_gc = SA[X[1], X[2], X[3], vpar]
                xv_rec = TP._gc_to_full(
                    state_gc, E_field, B_field_flat, q, m, μ, 0.0, phase
                )
                @test norm(xv[SA[1, 2, 3]] - xv_rec[SA[1, 2, 3]]) < 1.0e-10
                @test norm(xv[SA[4, 5, 6]] - xv_rec[SA[4, 5, 6]]) < 1.0e-8
            end
        end
    end

    # 5. Forced pure-FO hybrid must reproduce standalone Boris bit-for-bit.
    # (threshold = 0 forces FO for the whole trajectory since ε ≥ 0 always.)
    @testset "full orbit matches Boris" begin
        v0 = SA[5.0e4, 0.0, 1.0e5]
        u0 = vcat(SA[0.0, 0.0, 0.0], v0)
        T_gyro = gyroperiod(0.01)
        tspan = (0.0, 30 * T_gyro)
        p = (q2m, m, E_field, B_field_flat, TP.ZeroField())

        alg_fo = AdaptiveHybrid(;
            threshold = 0.0, dtmax = 10 * T_gyro, dtmin = 1.0e-6 * T_gyro,
            check_interval = 100
        )
        sol_h = solve(TraceProblem(u0, tspan, p), alg_fo).u[1]

        dt_fo = 2π * alg_fo.safety_fo / (abs(q2m) * 0.01)
        sol_b = solve(TraceProblem(u0, tspan, p), Boris(); dt = dt_fo)

        @test length(sol_h.t) == length(sol_b.t)
        for i in 1:length(sol_b.t)
            uh = sol_h(sol_b.t[i])
            ub = sol_b.u[i]
            @test isapprox(norm(uh[SA[1, 2, 3]] - ub[SA[1, 2, 3]]), 0.0; atol = 1.0e-8)
            @test isapprox(norm(uh[SA[4, 5, 6]] - ub[SA[4, 5, 6]]), 0.0; atol = 1.0e-5)
        end
    end

    # 6. Hysteresis thresholds: backward-compatible `threshold` sets both bounds,
    # and a wider buffer (α > β) reduces GC <-> FO chattering.
    @testset "hysteresis thresholds" begin
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0).threshold_gc_to_fo == 0.1
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0).threshold_fo_to_gc == 0.1

        # α must be >= β (buffer region requires α ≥ β).
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold_gc_to_fo = 0.05, threshold_fo_to_gc = 0.2, dtmax = 1.0
        )

        # Distinct bounds are stored verbatim.
        alg_buf = AdaptiveHybrid(;
            threshold_gc_to_fo = 0.2, threshold_fo_to_gc = 0.05, dtmax = 1.0
        )
        @test alg_buf.threshold_gc_to_fo == 0.2
        @test alg_buf.threshold_fo_to_gc == 0.05
    end

    # 7. Buffer region suppresses mode chatter: a particle whose adiabaticity
    # parameter oscillates within [β, α) must keep a single, stable mode when a
    # wide hysteresis band is used, whereas a zero-width band (α == β) permits
    # switching at every crossing. We exercise this with the bottle field and
    # verify the run still completes successfully with a wide buffer.
    @testset "wide buffer" begin
        u0 = SA[0.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        T_gyro = gyroperiod(bottle_B0)
        tspan = (0.0, 30 * T_gyro)
        p = (q2m, m, E_field, B_field_bottle, TP.ZeroField())

        # Wide hysteresis band: α = 0.2, β = 0.05.
        alg_buf = AdaptiveHybrid(;
            threshold_gc_to_fo = 0.2, threshold_fo_to_gc = 0.05,
            dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro, check_interval = 100,
        )
        sol_buf = solve(TraceProblem(u0, tspan, p), alg_buf).u[1]
        @test sol_buf.retcode == TP.ReturnCode.Success
        @test length(sol_buf.t) > 100
    end

    # 8. Saving at requested times, for a case that stays in the GC mode.
    @testset "saveat" begin
        u0 = vcat(SA[0.0, 0.0, 0.0], SA[1.0e4, 0.0, 1.0e3])
        tspan = (0.0, 1.0e-4)
        p = (q2m, m, E_field, B_field_flat, TP.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)
        prob = TraceProblem(u0, tspan, p)

        sol_all = solve(prob, alg).u[1]

        # Requesting intermediate output must not change the integration, only
        # where the state is reported.
        ts = collect(0.1e-4:0.1e-4:0.9e-4)
        sol_at = solve(prob, alg; saveat = ts).u[1]

        @test sol_at.t ≈ vcat(0.0, ts, 1.0e-4)
        @test sol_at.u[end] ≈ sol_all.u[end]

        # A shared time must give the same state whichever grid it arrives in,
        # which pins down that each one is interpolated inside its own step.
        ts_fine = sort!(vcat(ts, collect(0.05e-4:0.1e-4:0.95e-4)))
        sol_fine = solve(prob, alg; saveat = ts_fine).u[1]

        for t in ts
            j = findfirst(isequal(t), sol_fine.t)
            k = findfirst(isequal(t), sol_at.t)
            @test j !== nothing
            @test k !== nothing
            @test sol_fine.u[j] ≈ sol_at.u[k]
        end

        # A plain interval is accepted too: its interior is used, since the ends
        # of the span are reported separately.
        sol_interval = solve(prob, alg; saveat = 0.25e-4).u[1]
        @test sol_interval.t ≈ collect(0.0:0.25e-4:1.0e-4)
    end

    # 9. Constructor and parameter validation
    @testset "validation" begin
        p_dummy = (q2m, m, E_field, B_field_flat, TP.ZeroField())
        u0_6 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        @test_throws ArgumentError AdaptiveHybrid(; threshold = 0.1, dtmax = -1.0)
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, dtmin = -0.1
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, dtmin = 2.0
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, safety_fo = 0.0
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, maxiters = 0
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold_fo_to_gc = -0.1, threshold_gc_to_fo = 0.1, dtmax = 1.0
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, abstol = -1.0e-6
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, reltol = 0.0
        )

        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0)
        @test_throws ArgumentError solve(
            TraceProblem(SA[0.0, 0.0, 0.0], (0.0, 1.0), p_dummy), alg
        )
        @test_throws ArgumentError solve(TraceProblem(u0_6, (1.0, 0.0), p_dummy), alg)
    end

    # 10. Early termination with isoutside and maxiters
    @testset "early termination" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        tspan = (0.0, 1.0e-3)
        p = (q2m, m, E_field, B_field_flat, TP.ZeroField())
        prob = TraceProblem(u0, tspan, p)
        isoutside = (u, p, t) -> u[3] > 1.0e-2

        # Termination via isoutside in GC mode
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-5)
        sol_out = solve(prob, alg; isoutside).u[1]
        @test sol_out.retcode == TP.ReturnCode.Terminated
        @test sol_out.t[end] < tspan[2]

        # Termination via isoutside in FO mode (threshold = 0.0 forces FO)
        alg_fo = AdaptiveHybrid(; threshold = 0.0, dtmax = 1.0e-5)
        sol_out_fo = solve(prob, alg_fo; isoutside).u[1]
        @test sol_out_fo.retcode == TP.ReturnCode.Terminated
        @test sol_out_fo.t[end] < tspan[2]

        # Termination via maxiters
        alg_max = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6, maxiters = 3)
        sol_max = solve(prob, alg_max).u[1]
        @test sol_max.retcode == TP.ReturnCode.MaxIters
    end

    # 11. Minimal storage when save_everystep = false
    @testset "minimal storage" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        tspan = (0.0, 1.0e-4)
        p = (q2m, m, E_field, B_field_flat, TP.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)
        prob = TraceProblem(u0, tspan, p)

        sol_ends = solve(prob, alg; save_everystep = false).u[1]
        @test length(sol_ends.t) == 2
        @test sol_ends.t[1] == tspan[1]
        @test sol_ends.t[end] == tspan[2]
    end

    # 12. Monotonicity and transition preservation in hybrid switching
    @testset "diagnostic timestamps" begin
        u0 = SA[0.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        T_gyro = gyroperiod(bottle_B0)
        tspan = (0.0, 10 * T_gyro)
        p = (q2m, m, E_field, B_field_bottle, TP.ZeroField())
        alg = AdaptiveHybrid(;
            threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            check_interval = 10,
        )
        sol = solve(TraceProblem(u0, tspan, p), alg).u[1]
        @test sol.retcode == TP.ReturnCode.Success

        # Diagnostic timestamps must be strictly monotonic (no duplicates)
        adia = sol.stats.adiabaticity
        @test length(adia.t) > 1
        @test all(diff(adia.t) .> 0)
    end

    # 13. The full-orbit step follows the local gyroperiod instead of staying at
    # the value it had when the solver entered the full-orbit mode.
    @testset "adaptive full-orbit step" begin
        u0 = SA[0.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        T_gyro = gyroperiod(bottle_B0)
        p = (q2m, m, E_field, B_field_bottle, TP.ZeroField())

        # threshold = 0 forces the full-orbit mode for the whole trajectory.
        alg = AdaptiveHybrid(;
            threshold = 0.0, dtmax = T_gyro, dtmin = 1.0e-6 * T_gyro,
            check_interval = 10,
        )
        sol = solve(TraceProblem(u0, (0.0, 10 * T_gyro), p), alg).u[1]
        @test sol.retcode == TP.ReturnCode.Success

        steps = diff(sol.t)
        @test length(steps) > 10
        # The bottle field spans a factor of two in |B| over the bounce, so the
        # step size must change by more than a few percent. The last step is the
        # truncated remainder of the time span and is excluded.
        @test maximum(steps[1:(end - 1)]) > 1.05 * minimum(steps[1:(end - 1)])
    end

    # 14. The adiabaticity check is bounded in time, not only in accepted steps,
    # so a long guiding-center step cannot skip far into a non-adiabatic region.
    @testset "check time bound" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        T_gyro = gyroperiod(0.01)
        p = (q2m, m, E_field, B_field_flat, TP.ZeroField())
        prob = TraceProblem(u0, (0.0, 1.0e-2), p)

        # A uniform field keeps the guiding-center steps growing geometrically,
        # so the step count never reaches `check_interval` in this span.
        base = (;
            threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            check_interval = 1000,
        )
        sol_bound = solve(prob, AdaptiveHybrid(; base..., check_dt = T_gyro)).u[1]
        sol_free = solve(prob, AdaptiveHybrid(; base..., check_dt = 1.0e9)).u[1]

        @test sol_bound.retcode == TP.ReturnCode.Success
        @test sol_free.retcode == TP.ReturnCode.Success
        # With a time bound the solver checks repeatedly; without one only the
        # initial diagnostic is recorded.
        @test length(sol_bound.stats.adiabaticity.t) > 2
        @test length(sol_free.stats.adiabaticity.t) == 1
    end

    # 15. The switch correction refines the guiding-center map with fixed-point
    # iterations, and `switch_correction = 0` recovers the zeroth-order map.
    @testset "switch correction" begin
        T_gyro = gyroperiod(bottle_B0)
        # A state off the magnetic axis, where the field is curved.
        u = SA[5.0, 0.0, 40.0, 4.0e4, 1.0e4, 9.0e4]
        t = 0.0

        X0, vpar0, μ0, phase0 = TP._gc_state_at_switch(
            u, E_field, B_field_bottle, q, m, t, 0
        )
        X_legacy, vpar_legacy, _, _, μ_legacy, phase_legacy = TP._get_gc_parameters(
            u, E_field, B_field_bottle, q, m, t
        )
        @test X0 == X_legacy
        @test vpar0 == vpar_legacy
        @test μ0 == μ_legacy

        # Residual of the implicit map X = x - ρ(X): the correction reduces it
        # by orders of magnitude for a curved field.
        residual(X, vpar, μ, phase) = begin
            x_rec = TP._gc_to_full(
                SVector(X[1], X[2], X[3], vpar), E_field, B_field_bottle, q, m,
                μ, t, phase
            )
            norm(TP.get_x(x_rec) - u[SA[1, 2, 3]])
        end
        X2, vpar2, μ2, phase2 = TP._gc_state_at_switch(
            u, E_field, B_field_bottle, q, m, t, 2
        )
        @test residual(X2, vpar2, μ2, phase2) < 1.0e-2 * residual(X0, vpar0, μ0, phase0)

        # The corrected map is the exact inverse of the reconstruction used when
        # switching back to the full orbit: GC -> FO -> GC returns the state.
        state_gc = SVector(X2[1], X2[2], X2[3], vpar2)
        xv_back = TP._gc_to_full(
            state_gc, E_field, B_field_bottle, q, m, μ2, t, phase2
        )
        X_back, vpar_back, μ_back, phase_back = TP._gc_state_at_switch(
            xv_back, E_field, B_field_bottle, q, m, t, 2
        )
        @test X_back ≈ X2 atol = 1.0e-5
        @test vpar_back ≈ vpar2 rtol = 1.0e-6
        @test μ_back ≈ μ2 rtol = 1.0e-6

        # `switch_correction = 0` reproduces the legacy conversion exactly.
        alg = AdaptiveHybrid(;
            threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            switch_correction = 0,
        )
        sol = solve(
            TraceProblem(
                u, (0.0, 5 * T_gyro),
                (q2m, m, E_field, B_field_bottle, TP.ZeroField())
            ),
            alg
        ).u[1]
        @test sol.retcode == TP.ReturnCode.Success
    end

    # 16. Validation of the new keywords.
    @testset "new keyword validation" begin
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 2.0).check_dt == 20.0
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 2.0).switch_correction == 2
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, check_dt = 0.0
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, switch_correction = -1
        )
    end
end

end # module test_hybrid
