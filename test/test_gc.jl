module test_gc

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using TestParticle
import TestParticle as TP
using StaticArrays
using Test
using LinearAlgebra
using SciMLBase
using OrdinaryDiffEq
import Magnetostatics as MS
using ..test_common: uniform_Ex, curved_B

const Ek = 5.0e7 # [eV], for the dipole test case

"Earth's dipole field, the workhorse for the native GC solvers."
const dipole = MS.Dipole(TP.BMoment_Earth)

dipole_field(r) = dipole(r)

"Constant electric field with the time-dependent signature."
E_const(x, t) = SA[1.0e-9, 0.0, 0.0]

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

        # The shift of the guiding center is the Larmor radius, which goes as
        # 1/B, so the two cases have to differ by B_7 / B_6 = 2.
        norm_shift_7 = norm(gc_7[1:3])
        norm_shift_6 = norm(gc_6[1:3])

        @test isapprox(norm_shift_6 / norm_shift_7, 2.0, rtol = 0.01)
    end

    @testset "GC <-> Full Conversion" begin
        # Setup simple uniform field
        B0 = 1.0 # T
        B_field(r) = SA[0, 0, B0]
        E_field(r) = SA[0, 0, 0]
        param = prepare(E_field, B_field; species = Proton)

        @testset "Reversibility GC -> Full -> GC" begin
            R = SA[1.0, 0.0, 0.0]
            vpar = 1.0e5
            gc_state = [R..., vpar]
            μ = 1.0e-15
            phase = π / 3

            xu = gc_to_full(gc_state, param, μ, phase)
            @test length(xu) == 6

            gc_new, μ_new = full_to_gc(xu, param)

            @test gc_new[1:3] ≈ R rtol = 1.0e-5
            @test gc_new[4] ≈ vpar rtol = 1.0e-5
            @test μ_new ≈ μ rtol = 1.0e-5
        end

        @testset "With E-field (ExB drift)" begin
            # B in Z, E in Y -> ExB in X
            B0 = 1.0
            E0 = 1.0e3
            B_field2(r) = SA[0, 0, B0]
            E_field2(r) = SA[0, E0, 0] # E in y
            param2 = prepare(E_field2, B_field2; species = Proton)

            # v_E = (E x B) / B² = (E0/B0, 0, 0) = (1000, 0, 0)
            gc_state = [SA[0.0, 0.0, 0.0]..., 0.0]
            μ = 0.0 # No gyromotion, only drift

            xu = gc_to_full(gc_state, param2, μ)
            v = xu[4:6]
            @test v[1] ≈ 1000.0 atol = 1.0e-5
            @test v[2] ≈ 0.0 atol = 1.0e-5
            @test v[3] ≈ 0.0 atol = 1.0e-5

            gc_new, μ_new = full_to_gc(xu, param2)
            @test gc_new[1:3] ≈ SA[0.0, 0.0, 0.0] atol = 1.0e-5
            @test gc_new[4] ≈ 0.0 atol = 1.0e-5
            @test μ_new ≈ 0.0 atol = 1.0e-10
        end

        @testset "Phase check" begin
            B_field3(r) = SA[0, 0, 1.0]
            param3 = prepare(E_field, B_field3; species = Proton)

            gc_state = [SA[0.0, 0.0, 0.0]..., 0.0]
            μ = 1.0e-15

            # For B in Z the perpendicular plane is XY, so phase 0 points along
            # +x and phase π/2 along +y.
            v0 = gc_to_full(gc_state, param3, μ, 0.0)[4:6]
            v90 = gc_to_full(gc_state, param3, μ, π / 2)[4:6]

            @test dot(normalize(v0), normalize(v90)) ≈ 0.0 atol = 1.0e-5
            @test isapprox(normalize(v0), SA[1.0, 0.0, 0.0], atol = 1.0e-5)
            @test isapprox(normalize(v90), SA[0.0, 1.0, 0.0], atol = 1.0e-5)
        end
    end

    @testset "Native Solvers" begin
        # Setup simple dipole field
        m = TP.mᵢ
        q = TP.qᵢ
        Rₑ = TP.Rₑ

        # Initial condition
        v₀ = sph2cart(energy2velocity(Ek; q, m), π / 4, 0.0)
        r₀ = sph2cart(2.5 * Rₑ, π / 2, 0.0)
        stateinit = [r₀..., v₀...]
        tspan = (0.0, 1.0)

        stateinit_gc, param_gc = prepare_gc(
            stateinit, TP.ZeroField(), dipole_field;
            species = Proton
        )

        prob = TraceGCProblem(stateinit_gc, tspan, param_gc)

        @testset "Fixed RK4" begin
            dt = 1.0e-4
            sol = solve(prob; dt, alg = :rk4)
            @test length(sol.u) == 1
            @test length(sol.u[1].t) == 10001
            @test sol.u[1].retcode == ReturnCode.MaxIters

            # Test save_everystep=false
            sol_no_save = solve(prob; dt, alg = :rk4, save_everystep = false)
            @test length(sol_no_save.u[1].t) == 2 # start and end

            # Test early exit to cover resize!
            sol_early = solve(
                prob; dt, alg = :rk4,
                isoutside = (xv, p, t) -> t > 0.5
            )
            @test length(sol_early.u[1].t) == 5001
        end

        @testset "Adaptive RK45" begin
            # Default tolerances
            sol_def = solve(prob; dt = 1.0e-4, alg = :rk45)
            @test length(sol_def.u[1].t) == 61
            @test sol_def.u[1].retcode == ReturnCode.Success

            # Tight tolerances
            sol_tight = solve(
                prob;
                dt = 1.0e-4, alg = :rk45, abstol = 1.0e-8, reltol = 1.0e-8
            )
            @test length(sol_tight.u[1].t) > length(sol_def.u[1].t)

            # Accuracy check
            sol_rk4 = solve(prob; dt = 1.0e-4, alg = :rk4)
            diff = norm(sol_tight.u[1].u[end] - sol_rk4.u[1].u[end])
            @test diff < 10.0
        end

        @testset "Comparison with DiffEq" begin
            # Solve using DiffEq (Vern9, high accuracy)
            prob_diffeq = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
            sol_diffeq = solve(prob_diffeq, Vern9(); reltol = 1.0e-8, abstol = 1.0e-8)
            u_diffeq = sol_diffeq.u[end]

            # Solve using Native RK4
            dt = 1.0e-4
            sol_native = solve(prob; dt, alg = :rk4)
            u_native = sol_native.u[1].u[end]

            # Position difference
            @test norm(u_diffeq[1:3] - u_native[1:3]) / norm(u_diffeq[1:3]) < 1.0e-3

            # Parallel velocity difference
            @test abs(u_diffeq[4] - u_native[4]) / abs(u_diffeq[4]) < 1.0e-3
        end

        @testset "Automatic Initial dt" begin
            # Test that we can call solve without dt for adaptive method
            sol_auto = solve(prob; alg = :rk45)
            @test sol_auto.u[1].retcode == ReturnCode.Success
            @test length(sol_auto.u[1].t) == 58
        end

        @testset "Ensemble" begin
            trajectories = 10

            # Serial
            sol_serial = solve(prob; trajectories, dt = 1.0e-4, alg = :rk45)
            @test length(sol_serial.u) == trajectories
            @test all(s -> s.retcode == ReturnCode.Success, sol_serial.u)

            # Threads
            sol_threads = solve(
                prob, EnsembleThreads(); trajectories, dt = 1.0e-4, alg = :rk45
            )
            @test length(sol_threads.u) == trajectories
            @test all(s -> s.retcode == ReturnCode.Success, sol_threads.u)
        end

        @testset "GC Extra Saving" begin
            trajectories = 2
            tspan_saving = (0.0, 1.0e-4)
            r0 = [1.5 * TP.Rₑ, 0.0, 0.0]
            v0 = [0.0, 1.0e5, 1.0e5] # Arbitrary
            state_gc_0, param_gc_ready = prepare_gc(
                vcat(r0, v0), TP.ZeroField(), dipole_field; species = Proton
            )
            prob_saving = TraceGCProblem(state_gc_0, tspan_saving, param_gc_ready)

            # Solve with saving enabled
            sol = solve(
                prob_saving; trajectories, dt = 1.0e-5, save_fields = true, save_work = true
            )

            for s in sol.u
                E, B = get_fields(s)
                W = get_work(s)

                @test length(E) == length(s.t)
                @test length(B) == length(s.t)
                @test length(W) == length(s.t)

                # Check types and basic values
                @test E[1] isa SVector{3}
                @test B[1] isa SVector{3}
                @test W[1] isa SVector{4}

                # Dipole field should have B magnitude non-zero
                @test norm(B[1]) > 0.0
                # Zero E-field
                @test norm(E[1]) == 0.0
            end
        end

        @testset "saveat" begin
            dt = 1.0e-4
            sol_all = solve(prob; dt, alg = :rk4).u[1]

            # Asking for intermediate output must not change the integration,
            # only where the state is reported. `save_start` and `save_end` keep
            # adding the two ends of the span around the requested times.
            ts = collect(0.1:0.1:0.9)
            sol_at = solve(prob; dt, alg = :rk4, saveat = ts).u[1]

            @test sol_at.t ≈ vcat(0.0, ts, 1.0)
            for (k, t) in enumerate(sol_at.t)
                @test sol_at.u[k] ≈ sol_all(t)
            end

            # The adaptive solver takes the same times.
            sol_at45 = solve(prob; dt, alg = :rk45, saveat = ts).u[1]
            @test sol_at45.t ≈ vcat(0.0, ts, 1.0)

            # An interval is accepted as well, and only the interior is taken
            # from it since the ends are reported separately.
            sol_interval = solve(prob; dt, alg = :rk4, saveat = 0.25).u[1]
            @test sol_interval.t ≈ collect(0.0:0.25:1.0)

            # The field and work columns are appended on this path too.
            sol_fields = solve(
                prob; dt, alg = :rk4, saveat = ts, save_fields = true
            ).u[1]
            @test length(sol_fields.u[1]) == 10
            @test get_fields(sol_fields) isa Tuple
        end
    end
end

@testset "GC drifts" begin
    stateinit = [1.0, 0.0, 0.0, 0.0, 1.0, 0.1]
    tspan = (0, 10)

    param = prepare(uniform_Ex, curved_B, species = Proton)
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Vern9())

    @views begin
        x, y, z = sol[1, :], sol[2, :], sol[3, :]
        vx, vy, vz = sol[4, :], sol[5, :], sol[6, :]
    end
    b = [curved_B(SA[xi, yi, zi]) for (xi, yi, zi) in zip(x, y, z)]
    bx = [bi[1] for bi in b]
    by = [bi[2] for bi in b]
    bz = [bi[3] for bi in b]
    X = get_gc(x, y, z, vx, vy, vz, bx, by, bz, param[1])
    @test sum(X[end]) ≈ 1.5380318687643348 rtol = 1.0e-6

    stateinit_gc,
        param_gc = prepare_gc(
        stateinit, uniform_Ex, curved_B,
        species = Proton
    )

    prob_gc = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
    sol_gc = solve(prob_gc, Vern9())

    # analytical drifts
    gc = param |> get_gc_func
    gc_x0 = gc(stateinit) |> Vector # needs mutation
    prob_gc_analytic = ODEProblem(trace_gc_drifts!, gc_x0, tspan, (param..., sol))
    sol_gc_analytic = solve(prob_gc_analytic, Vern9(); save_idxs = [1, 2, 3])
    @test sol_gc[1, end] ≈ 0.9896197928850492
    @test sol_gc_analytic[1, end] ≈ 0.9906923500002904 rtol = 1.0e-5

    # Test get_E_parameters with constant E field (not used currently)
    x_test, t_test = SA[1.0, 0.0, 0.0], 0.0
    E_expected = E_const(x_test, t_test)
    E, JE = TP.get_E_parameters(x_test, t_test, E_const)
    @test E == E_expected && JE == zeros(3, 3)
end

end # module test_gc
