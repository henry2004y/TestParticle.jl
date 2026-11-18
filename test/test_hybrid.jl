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
using OrdinaryDiffEq
using OrdinaryDiffEqVerner
using Test
using ..test_common: bottle_B, bottle_B0, sheared_B

const m = TP.mᵢ
const q = TP.qᵢ
const q2m = q / m

const E_field = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

"Uniform 0.01 T field along +z, flat enough to stay adiabatic."
B_flat(x, t) = SA[0.0, 0.0, 0.01]

"Harris-type current sheet with weak normal guide component."
function harris_B(x, t)
    B0 = 2.0e-8
    Bn = 2.0e-9
    L = 2.0e4
    return SA[B0 * tanh(x[3] / L), 0.0, Bn]
end

"Thin reversing sheet with weak normal field so X_gc teleports far into uniform lobes."
near_null_sheet_B(x, t) = SA[1.0e-5 * tanh(x[3] / 1.0e3), 0.0, 1.0e-11]

"Thin reversing sheet with an exact magnetic null at z = 0."
exact_null_sheet_B(x, t) = SA[1.0e-5 * tanh(x[3] / 1.0e3), 0.0, 0.0]

"Adiabatic magnetic bottle (α = 1e-8) for testing step preservation."
function adia_bottle_B(x, t)
    B0 = 1.0e-4
    α = 1.0e-8
    return SA[-B0 * α * x[1] * x[3], -B0 * α * x[2] * x[3], B0 * (1 + α * x[3]^2)]
end

"Harris current sheet from benchmark suite for overshoot/reflection checks."
sheet_bench_B(x, t) = SA[1.0e-7 * tanh(x[3] / 2.0e4), 0.0, 2.0e-8]

"Float32-preserving sheared magnetic field."
function sheared_B32(x, t)
    T = eltype(x)
    return SA[T(1.0e-4), T(5.0e-4) * x[3], T(0)]
end

"Float32-preserving magnetic bottle."
function bottle_B32(x, t)
    T = eltype(x)
    B0 = T(1.0e-6)
    L = T(100.0)
    k = one(T) + (x[3] / L)^2
    dk = T(2) * x[3] / L^2
    return SA[-T(0.5) * B0 * x[1] * dk, -T(0.5) * B0 * x[2] * dk, B0 * k]
end

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
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
        alg = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6)

        sols = solve(TraceProblem(u0, tspan, p), alg)
        @test sols.u[1].retcode == TP.ReturnCode.Success
    end

    # 2. Non-adiabatic Case (Highly Curved Field)
    @testset "sheared field" begin
        u0 = vcat(SA[0.0, 0.0, 0.0], SA[1.0e4, 0.0, 1.0e3])
        tspan = (0.0, 1.0e-5)
        p = (; q2m, m, E = E_field, B = B_field_sheared, F = TP.ZeroField())
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
        p = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())

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
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())

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
        p = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())

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
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
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
        p_dummy = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
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
        @test_throws ArgumentError solve(TraceProblem(u0_6, (1.0, 1.0), p_dummy), alg)
    end

    # 10. Early termination with isoutside and maxiters
    @testset "early termination" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        tspan = (0.0, 1.0e-3)
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
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
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
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
        p = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())
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
        p = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())

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

    # 14. The step never crosses the next check deadline, so a long
    # guiding-center step cannot pass through a non-adiabatic region unexamined.
    @testset "check deadline" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        T_gyro = gyroperiod(0.01)
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
        check_dt = 5 * T_gyro
        alg = AdaptiveHybrid(;
            threshold = 0.1, dtmax = 10 * T_gyro, dtmin = 1.0e-6 * T_gyro,
            check_dt, check_interval = 1_000_000,
        )
        sol = solve(TraceProblem(u0, (0.0, 40 * T_gyro), p), alg).u[1]
        @test sol.retcode == TP.ReturnCode.Success

        gaps = diff(sol.stats.adiabaticity.t)
        push!(gaps, sol.t[end] - sol.stats.adiabaticity.t[end])
        @test length(gaps) > 2
        @test maximum(gaps) <= check_dt * (1 + 1.0e-9)
    end

    # 14a. Landing on a deadline repeatedly must still trigger a check: an
    # elapsed-time comparison can round just below `check_dt` and skip it.
    @testset "repeated deadline landings" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        T_gyro = gyroperiod(0.01)
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
        for (threshold, check_dt, dtmin, dtmax, tspan) in (
                (0.1, 0.5 * T_gyro, 1.0e-6 * T_gyro, 10 * T_gyro, 40 * T_gyro),
                (0.0, 0.5 * T_gyro, 1.0e-6 * T_gyro, 10 * T_gyro, 40 * T_gyro),
                (0.1, 0.1, 1.0e-9, 1.0, 0.5),
            )
            alg = AdaptiveHybrid(;
                threshold, dtmax, dtmin, check_dt, check_interval = 1_000_000,
            )
            sol = solve(TraceProblem(u0, (0.0, tspan), p), alg).u[1]
            @test sol.retcode == TP.ReturnCode.Success

            ad = sol.stats.adiabaticity
            gaps = diff(ad.t)
            push!(gaps, sol.t[end] - ad.t[end])
            @test length(gaps) > 2
            @test maximum(gaps) <= check_dt * (1 + 1.0e-9)
        end
    end

    # 14b. A positive deadline remainder is clipped even when it is below
    # `dtmin`, which is an allowed configuration (`check_dt <= dtmin`).
    @testset "deadline below dtmin" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        T_gyro = gyroperiod(0.01)
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
        check_dt = 0.01 * T_gyro
        alg = AdaptiveHybrid(;
            threshold = 0.1, dtmax = 10 * T_gyro, dtmin = 0.1 * T_gyro,
            check_dt, check_interval = 1_000_000,
        )
        sol = solve(TraceProblem(u0, (0.0, 2 * T_gyro), p), alg).u[1]
        @test sol.retcode == TP.ReturnCode.Success

        ad = sol.stats.adiabaticity
        gaps = diff(ad.t)
        push!(gaps, sol.t[end] - ad.t[end])
        @test length(gaps) > 20
        @test maximum(gaps) <= check_dt * (1 + 1.0e-9)
    end

    # 15. The adiabaticity check is bounded in time, not only in accepted steps,
    # so a long guiding-center step cannot skip far into a non-adiabatic region.
    @testset "check time bound" begin
        u0 = SA[0.0, 0.0, 0.0, 1.0e4, 0.0, 1.0e3]
        T_gyro = gyroperiod(0.01)
        p = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
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
                (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())
            ),
            alg
        ).u[1]
        @test sol.retcode == TP.ReturnCode.Success
    end

    # 16. The adiabaticity check reuses the field derivatives that the step needs
    # anyway, so both must agree with the standalone evaluation.
    @testset "cached field properties" begin
        u = SA[5.0, 0.0, 40.0, 4.0e4, 1.0e4, 9.0e4]
        X = u[SA[1, 2, 3]]
        μ = 1.0e-15
        mp = TP.get_magnetic_properties(X, 0.0, B_field_bottle)

        @test all(
            isapprox.(
                TP._adiabaticity_components(mp, q, m, μ),
                TP.adiabaticity_components(X, B_field_bottle, q, m, μ, 0.0)
            )
        )

        # The paths compared below compile separately, and `muladd` is free to fuse
        # differently in each, so they agree only up to rounding on some CPUs.

        # The lean three-element version used by the check drops the gyrophase but
        # returns the same guiding center, parallel velocity, and magnetic moment.
        X_check, vpar_check, μ_check = TP._gc_adia_parameters(
            u, E_field, B_field_bottle, q, m, 0.0
        )
        X_full, vpar_full, _, _, μ_full, _ = TP._get_gc_parameters(
            u, E_field, B_field_bottle, q, m, 0.0
        )
        @test X_check ≈ X_full
        @test vpar_check ≈ vpar_full
        @test μ_check ≈ μ_full

        # Custom `gc_alg` can be passed to `AdaptiveHybrid`.
        u0_flat = vcat(SA[0.0, 0.0, 0.0], SA[1.0e4, 0.0, 1.0e3])
        p_flat = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())
        alg_tsit = AdaptiveHybrid(; threshold = 0.1, dtmax = 1.0e-6, gc_alg = Tsit5())
        sol_tsit = solve(TraceProblem(u0_flat, (0.0, 1.0e-4), p_flat), alg_tsit).u[1]
        @test sol_tsit.retcode == TP.ReturnCode.Success
    end

    # 17. A mode switch clips the step to the target mode's stability limit
    # instead of resetting it, so the step does not jump at the switch.
    @testset "step carried through a switch" begin
        u0 = SA[10.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        T_gyro = gyroperiod(bottle_B0)
        tspan = (0.0, 30 * T_gyro)
        p = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())
        alg = AdaptiveHybrid(;
            threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            check_interval = 10,
        )
        sol = solve(TraceProblem(u0, tspan, p), alg).u[1]
        @test sol.retcode == TP.ReturnCode.Success

        adia = sol.stats.adiabaticity
        steps = diff(sol.t)
        n_to_gc, n_to_fo = 0, 0
        for i in 2:length(adia.mode)
            adia.mode[i] === adia.mode[i - 1] && continue
            j = findfirst(t -> t ≈ adia.t[i], sol.t)
            (j === nothing || j < 2 || j >= length(steps)) && continue
            if adia.mode[i] === :GC
                # The full-orbit step is carried over, capped at half a
                # gyroperiod, so it never jumps above the cap.
                n_to_gc += 1
                @test steps[j] <= 1.2 * max(steps[j - 1], 0.5 * T_gyro)
            else
                # The guiding-center step is clipped to the full-orbit limit.
                n_to_fo += 1
                Bmag = norm(bottle_B(sol.u[j][SA[1, 2, 3]], sol.t[j]))
                dt_lim = min(2π * alg.safety_fo / (abs(q2m) * Bmag), alg.dtmax)
                @test steps[j] <= 1.2 * dt_lim
            end
        end
        @test n_to_gc > 0
        @test n_to_fo > 0
    end

    # 18. While the solver is confidently non-adiabatic, the full-orbit check
    # interval is stretched, so a run in that regime needs fewer checks.
    @testset "stretched check interval" begin
        u0 = SA[10.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        T_gyro = gyroperiod(bottle_B0)
        p = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())
        prob = TraceProblem(u0, (0.0, 30 * T_gyro), p)
        base = (;
            threshold = 0.05, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            check_interval = 10,
        )
        sol_flat = solve(prob, AdaptiveHybrid(; base..., check_margin_max = 1.0)).u[1]
        sol_long = solve(prob, AdaptiveHybrid(; base..., check_margin_max = 10.0)).u[1]

        @test sol_flat.retcode == TP.ReturnCode.Success
        @test sol_long.retcode == TP.ReturnCode.Success
        @test sol_long.stats.adiabaticity.t[2] > 2 * sol_flat.stats.adiabaticity.t[2]

        u0_sh = SA[0.0, 0.0, 0.0, 1.0e5, 0.0, 1.0e5]
        p_sh = (; q2m, m, E = E_field, B = B_field_sheared, F = TP.ZeroField())
        prob_sh = TraceProblem(u0_sh, (0.0, 1.0e-4), p_sh)
        base_sh = (;
            threshold = 0.1, dtmax = 1.0e-6, adiabaticity = :jacobian,
            check_interval = 10,
        )
        sol_sh_flat = solve(
            prob_sh, AdaptiveHybrid(; base_sh..., check_margin_max = 1.0)
        ).u[1]
        sol_sh_long = solve(
            prob_sh, AdaptiveHybrid(; base_sh..., check_margin_max = 10.0)
        ).u[1]
        @test sol_sh_flat.retcode == TP.ReturnCode.Success
        @test sol_sh_long.retcode == TP.ReturnCode.Success
        @test length(sol_sh_long.stats.adiabaticity.t) <
            length(sol_sh_flat.stats.adiabaticity.t)
    end

    # 19. Validation of the new keywords.
    @testset "new keyword validation" begin
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 2.0).check_dt == 200.0
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 2.0).switch_correction == 2
        @test AdaptiveHybrid(; threshold = 0.1, dtmax = 2.0).check_margin_max == 1.0
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, check_dt = 0.0
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, switch_correction = -1
        )
        @test_throws ArgumentError AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0, check_margin_max = 0.5
        )
    end

    @testset "field-aligned current sheet" begin
        L = 2.0e4
        r0 = SA[0.0, 0.0, -5 * L]
        B_init = harris_B(r0, 0.0)
        b_init = B_init / sqrt(sum(abs2, B_init))
        v0 = 4.0e5 * b_init
        u0 = vcat(r0, v0)
        p = (; q2m, m, E = E_field, B = TP.Field(harris_B), F = TP.ZeroField())
        prob = TraceProblem(u0, (0.0, 20.0), p)
        alg = AdaptiveHybrid(; threshold = 0.05, dtmax = 0.05, check_interval = 5)
        sol_hyb = solve(prob, alg).u[1]
        sol_boris = solve(prob, Boris(); dt = 0.05)

        @test sol_hyb.retcode == TP.ReturnCode.Success
        @test :GC in sol_hyb.stats.adiabaticity.mode
        @test :FO in sol_hyb.stats.adiabaticity.mode
        _, _, μ_init = TP._gc_adia_parameters(u0, E_field, p.B, q, m, 0.0)
        _, _, μ_boris = TP._gc_adia_parameters(
            sol_boris.u[end], E_field, p.B, q, m, sol_boris.t[end]
        )
        _, _, μ_hyb = TP._gc_adia_parameters(
            sol_hyb.u[end], E_field, p.B, q, m, sol_hyb.t[end]
        )
        @test μ_init < 1.0e-30
        @test μ_boris > 1.0e-13
        @test μ_hyb > 1.0e-13
        @test sign(sol_hyb.u[end][3]) == sign(sol_boris.u[end][3])
    end

    @testset "current-sheet null guard" begin
        # Near-null sheet where candidate X_gc teleports into uniform lobe (ε(X_gc) ≈ 0),
        # while local ε(r) is huge.
        u0_near = SA[0.0, 0.0, 1.0, 0.0, 1.0e5, 1.0e3]
        p_near = (;
            q2m, m, E = E_field, B = TP.Field(near_null_sheet_B), F = TP.ZeroField(),
        )
        alg_near = AdaptiveHybrid(;
            threshold = 0.1, dtmax = 0.01, check_interval = 2,
            adiabaticity = :both,
        )
        sol_near = solve(TraceProblem(u0_near, (0.0, 0.05), p_near), alg_near).u[1]
        @test sol_near.retcode == TP.ReturnCode.Success
        @test all(==(:FO), sol_near.stats.adiabaticity.mode)
        @test all(>(0.1), sol_near.stats.adiabaticity.components)

        # Exact null sheet crossing z = 0 and starting at z = 0 should produce no NaNs
        # and stay in :FO.
        u0_null = SA[0.0, 0.0, 1.0, 0.0, 1.0e5, -2.0e2]
        p_null = (;
            q2m, m, E = E_field, B = TP.Field(exact_null_sheet_B), F = TP.ZeroField(),
        )
        alg_null = AdaptiveHybrid(;
            threshold = 0.1, dtmax = 0.005, check_interval = 1,
            adiabaticity = :both,
        )
        sol_null = solve(TraceProblem(u0_null, (0.0, 0.02), p_null), alg_null).u[1]
        @test sol_null.retcode == TP.ReturnCode.Success
        @test all(==(:FO), sol_null.stats.adiabaticity.mode)
        @test !any(isnan, sol_null.stats.adiabaticity.components)
        @test !any(u -> any(isnan, u), sol_null.u)

        u0_zero = SA[0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0]
        sol_zero = solve(TraceProblem(u0_zero, (0.0, 0.02), p_null), alg_null).u[1]
        @test sol_zero.retcode == TP.ReturnCode.Success
        @test all(==(:FO), sol_zero.stats.adiabaticity.mode)
        @test !any(isnan, sol_zero.stats.adiabaticity.components)
        @test !any(u -> any(isnan, u), sol_zero.u)
    end

    @testset "saveat in FO and mixed modes" begin
        T_gyro = gyroperiod(bottle_B0)
        u0_fo = SA[10.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        p_bottle = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())
        ts_fo = collect(range(0.0, 5 * T_gyro; length = 6))
        alg_fo = AdaptiveHybrid(;
            threshold = 1.0e-5, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
        )
        sol_fo = solve(
            TraceProblem(u0_fo, (0.0, 5 * T_gyro), p_bottle), alg_fo; saveat = ts_fo
        ).u[1]
        @test sol_fo.retcode == TP.ReturnCode.Success
        @test all(==(:FO), sol_fo.stats.adiabaticity.mode)
        @test sol_fo.t ≈ ts_fo
        @test length(sol_fo.u) == length(ts_fo)

        ts_mix = collect(range(0.0, 30 * T_gyro; length = 16))
        alg_mix = AdaptiveHybrid(;
            threshold = 0.05, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            check_interval = 10,
        )
        sol_mix = solve(
            TraceProblem(u0_fo, (0.0, 30 * T_gyro), p_bottle), alg_mix; saveat = ts_mix
        ).u[1]
        @test sol_mix.retcode == TP.ReturnCode.Success
        @test :GC in sol_mix.stats.adiabaticity.mode
        @test :FO in sol_mix.stats.adiabaticity.mode
        @test sol_mix.t ≈ ts_mix
        @test length(sol_mix.u) == length(ts_mix)
    end

    @testset "backward integration" begin
        T_gyro = gyroperiod(bottle_B0)
        p_bottle = (; q2m, m, E = E_field, B = B_field_bottle, F = TP.ZeroField())
        p_flat = (; q2m, m, E = E_field, B = B_field_flat, F = TP.ZeroField())

        # Pure FO backward + reversibility
        u0_fo = SA[10.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        alg_fo = AdaptiveHybrid(;
            threshold = 1.0e-5, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            safety_fo = 0.005,
        )
        t_fo = 2 * T_gyro
        sol_fwd_fo = solve(TraceProblem(u0_fo, (0.0, t_fo), p_bottle), alg_fo).u[1]
        sol_bwd_fo = solve(
            TraceProblem(sol_fwd_fo.u[end], (t_fo, 0.0), p_bottle), alg_fo
        ).u[1]
        @test sol_bwd_fo.retcode == TP.ReturnCode.Success
        @test all(==(:FO), sol_bwd_fo.stats.adiabaticity.mode)
        @test sol_bwd_fo.t[1] == t_fo
        @test sol_bwd_fo.t[end] == 0.0
        @test issorted(sol_bwd_fo.t; rev = true)
        @test sol_bwd_fo.u[end][SA[1, 2, 3]] ≈ u0_fo[SA[1, 2, 3]] rtol = 1.0e-3
        @test sol_bwd_fo.u[end][SA[4, 5, 6]] ≈ u0_fo[SA[4, 5, 6]] rtol = 1.0e-3

        # Pure GC backward with saveat, small check_dt deadline clipping, and reversibility
        u0_gc = SA[1.0, 0.0, 0.0, 1.0e3, 0.0, 1.0e4]
        t_gc = 5.0e-4
        alg_gc = AdaptiveHybrid(;
            threshold = 0.1, dtmax = 1.0e-4, check_dt = 2.0e-5, check_interval = 100,
        )
        sol_fwd_gc = solve(TraceProblem(u0_gc, (0.0, t_gc), p_flat), alg_gc).u[1]
        ts_bwd_gc = collect(range(t_gc, 0.0; length = 6))
        sol_bwd_gc = solve(
            TraceProblem(sol_fwd_gc.u[end], (t_gc, 0.0), p_flat), alg_gc;
            saveat = ts_bwd_gc,
        ).u[1]
        @test sol_bwd_gc.retcode == TP.ReturnCode.Success
        @test all(==(:GC), sol_bwd_gc.stats.adiabaticity.mode)
        @test sol_bwd_gc.t ≈ ts_bwd_gc
        @test sol_bwd_gc.u[end][SA[1, 2, 3]] ≈ u0_gc[SA[1, 2, 3]] rtol = 1.0e-4
        @test sol_bwd_gc.u[end][SA[4, 5, 6]] ≈ u0_gc[SA[4, 5, 6]] rtol = 1.0e-4

        # Mixed FO/GC backward
        u0_mix = SA[10.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        t_mix = 30 * T_gyro
        alg_mix = AdaptiveHybrid(;
            threshold = 0.05, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro, check_interval = 10,
        )
        sol_bwd_mix = solve(TraceProblem(u0_mix, (t_mix, 0.0), p_bottle), alg_mix).u[1]
        @test sol_bwd_mix.retcode == TP.ReturnCode.Success
        @test sol_bwd_mix.t[1] == t_mix
        @test sol_bwd_mix.t[end] == 0.0
        @test issorted(sol_bwd_mix.t; rev = true)
        @test :GC in sol_bwd_mix.stats.adiabaticity.mode
        @test :FO in sol_bwd_mix.stats.adiabaticity.mode
    end

    @testset "Float32 preservation" begin
        # Electron in sheared field (ensures μ * m32 never underflows to zero)
        q_e32 = Float32(-abs(TP.qₑ))
        m_e32 = Float32(TP.mₑ)
        p_e32 = (;
            q2m = q_e32 / m_e32, m = m_e32, E = TP.ZeroField(),
            B = TP.Field(sheared_B32), F = TP.ZeroField(),
        )
        u0_e32 = Float32[0.0, 0.0, 0.0, 0.0, 1.0e5, 0.0]
        alg_e32 = AdaptiveHybrid(;
            threshold = 0.1f0, dtmax = 1.0f-9, adiabaticity = :jacobian,
        )
        sol_e32 = solve(TraceProblem(u0_e32, (0.0f0, 5.0f-9), p_e32), alg_e32).u[1]
        @test sol_e32.retcode == TP.ReturnCode.Success
        @test eltype(sol_e32.u[1]) == Float32
        @test eltype(sol_e32.t) == Float32
        @test eltype(sol_e32.stats.adiabaticity.t) == Float32
        @test eltype(sol_e32.stats.adiabaticity.components) == Float32
        @test !any(isnan, sol_e32.stats.adiabaticity.components)
        @test all(>(0.0f0), sol_e32.stats.adiabaticity.components)
        μ_e32 = Float32(0.5 * Float64(m_e32) * 1.0e10 / 1.0e-4)
        @test TP.adiabaticity_components(
            SA[0.0f0, 0.0f0, 0.0f0], sheared_B32, q_e32, m_e32, μ_e32, 0.0f0
        )[3] > 0.0f0

        # Proton in magnetic bottle
        q_p32 = Float32(q)
        m_p32 = Float32(m)
        T_p32 = Float32(gyroperiod(1.0e-6))
        p_p32 = (;
            q2m = q_p32 / m_p32, m = m_p32, E = TP.ZeroField(),
            B = TP.Field(bottle_B32), F = TP.ZeroField(),
        )
        u0_p32 = Float32[10.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
        alg_p32 = AdaptiveHybrid(;
            threshold = 0.05f0, dtmax = T_p32, dtmin = 1.0f-4 * T_p32,
        )
        sol_p32 = solve(TraceProblem(u0_p32, (0.0f0, 5 * T_p32), p_p32), alg_p32).u[1]
        @test sol_p32.retcode == TP.ReturnCode.Success
        @test eltype(sol_p32.u[1]) == Float32
        @test eltype(sol_p32.t) == Float32
    end

    # issue #667: check deadlines must preserve the GC integrator's proposed step
    # instead of crushing Vern7 back down to a remainder sliver.
    @testset "adiabatic step preservation" begin
        let
            T_gyro = gyroperiod(1.0e-4)
            u0 = SA[10.0, 0.0, 0.0, 5.0e4, 0.0, 1.0e5]
            tspan = (0.0, 20_000 * T_gyro)
            B_adia = TP.Field(adia_bottle_B)
            p = (; q2m, m, E = E_field, B = B_adia, F = TP.ZeroField())
            prob = TraceProblem(u0, tspan, p)

            u0_gc, p_gc = TP.prepare_gc(u0; B = B_adia, species = Proton)
            sol_gc = solve(
                TraceGCProblem(u0_gc, tspan, p_gc), Vern7();
                abstol = 1.0e-6, reltol = 1.0e-6,
            )

            alg_def = AdaptiveHybrid(;
                threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
                maxiters = 100_000,
            )
            sol_hyb = solve(prob, alg_def).u[1]
            @test sol_hyb.retcode == TP.ReturnCode.Success
            @test all(==(:GC), sol_hyb.stats.adiabaticity.mode)
            @test length(sol_hyb.t) <= ceil(Int, 1.15 * length(sol_gc.t))

            alg_clip = AdaptiveHybrid(;
                threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
                check_dt = 50 * T_gyro, maxiters = 100_000,
            )
            sol_clip = solve(prob, alg_clip).u[1]
            @test sol_clip.retcode == TP.ReturnCode.Success
            @test length(sol_clip.t) <= 430
        end
    end

    # issue #667: detect sheet entry early enough to reflect and capture the
    # non-adiabatic magnetic moment jump.
    @testset "current-sheet overshoot rollback" begin
        let
            B0 = 1.0e-7
            L = 2.0e4
            T_gyro = gyroperiod(B0)
            tspan = (0.0, 200 * T_gyro)
            v = 2.0e5
            B_sheet = TP.Field(sheet_bench_B)
            p = (; q2m, m, E = E_field, B = B_sheet, F = TP.ZeroField())
            dt_b = T_gyro / 20
            alg = AdaptiveHybrid(;
                threshold = 0.1, dtmax = dt_b, dtmin = 1.0e-4 * T_gyro,
                maxiters = 500_000, check_interval = 10,
            )

            u0 = vcat(SA[0.0, 0.0, -2 * L], SA[0.0, 0.95 * v, 0.3 * v])
            prob = TraceProblem(u0, tspan, p)
            sol_hyb = solve(prob, alg).u[1]
            sol_boris = solve(prob, Boris(); dt = dt_b)

            @test sol_hyb.retcode == TP.ReturnCode.Success
            idx_fo = findfirst(==(:FO), sol_hyb.stats.adiabaticity.mode)
            @test idx_fo !== nothing
            @test sol_hyb.stats.adiabaticity.components[idx_fo] <= 0.5
            @test sign(sol_hyb.u[end][3]) == sign(sol_boris.u[end][3]) == -1.0

            _, _, μ_boris = TP._gc_adia_parameters(
                sol_boris.u[end], E_field, B_sheet, q, m, tspan[2]
            )
            _, _, μ_hyb = TP._gc_adia_parameters(
                sol_hyb.u[end], E_field, B_sheet, q, m, tspan[2]
            )
            @test isapprox(μ_hyb, μ_boris; rtol = 0.2)

            # Starting deeper in the lobe exercises the GC -> FO crossing path.
            for s in (-10.5, -10.0, -9.95, -9.5)
                u0_lobe = vcat(SA[0.0, 0.0, s * L], SA[0.0, 0.95 * v, 0.3 * v])
                sol_lobe = solve(TraceProblem(u0_lobe, tspan, p), alg).u[1]
                @test sol_lobe.retcode == TP.ReturnCode.Success
                @test sol_lobe.stats.adiabaticity.mode[1] === :GC
                @test all(diff(sol_lobe.stats.adiabaticity.t) .> 0)
                idx_lobe_fo = findfirst(==(:FO), sol_lobe.stats.adiabaticity.mode)
                @test idx_lobe_fo !== nothing
                @test sol_lobe.stats.adiabaticity.components[idx_lobe_fo] <= 0.5
            end

            # Mode-switch info logging (`verbose = true`).
            sol_verb = @test_logs (:info, r"Initial mode") (:info, r"Switch") match_mode = :any begin
                solve(prob, alg; verbose = true).u[1]
            end
            @test sol_verb.retcode == TP.ReturnCode.Success

            # `isoutside` termination inside GC -> FO bisection.
            u0_s = SA[0.0, 0.0, -4.0 * L, 0.0, 0.95 * v, 0.3 * v]
            alg_iso = AdaptiveHybrid(;
                threshold = 0.05, dtmax = dt_b, dtmin = 1.0e-4 * T_gyro,
                check_dt = 50 * T_gyro, check_interval = 5,
            )
            sol_iso = solve(
                TraceProblem(u0_s, tspan, p), alg_iso;
                isoutside = (u, p, t) -> -1.3 * L < u[3] < -1.1 * L,
            ).u[1]
            @test sol_iso.retcode == TP.ReturnCode.Terminated

            # Direct switch at `t` when prior unchecked GC step sits in `[threshold, 2 * threshold]`.
            alg_prev = AdaptiveHybrid(;
                threshold = 0.15, dtmax = dt_b, dtmin = 1.0e-4 * T_gyro,
                check_dt = 50 * T_gyro, check_interval = 5,
            )
            sol_prev = solve(TraceProblem(u0_s, tspan, p), alg_prev).u[1]
            @test sol_prev.retcode == TP.ReturnCode.Success
            idx_prev = findfirst(==(:FO), sol_prev.stats.adiabaticity.mode)
            @test idx_prev !== nothing
            @test 0.15 <= sol_prev.stats.adiabaticity.components[idx_prev] <= 0.3

            # `maxiters` reached during GC -> FO bisection before reaching boundary.
            u0_close = SA[0.0, 0.0, -2.4 * L, 0.0, 0.95 * v, 0.3 * v]
            alg_mi = AdaptiveHybrid(;
                threshold = 0.1, dtmax = dt_b, dtmin = 1.0e-4 * T_gyro,
                check_dt = 50 * T_gyro, check_interval = 1, maxiters = 1,
            )
            sol_mi = solve(TraceProblem(u0_close, tspan, p), alg_mi).u[1]
            @test sol_mi.retcode == TP.ReturnCode.MaxIters
            @test all(==(:GC), sol_mi.stats.adiabaticity.mode)
        end
    end
end

end # module test_hybrid
