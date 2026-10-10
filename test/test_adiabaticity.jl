if !isdefined(Main, :test_common)
    include("test_common.jl")
end

module test_adiabaticity

using TestParticle
import TestParticle as TP
using StaticArrays
using LinearAlgebra
using Test
using ..test_common: bottle_B, bottle_B0, curved_B_t

const m = TP.mᵢ
const q = TP.qᵢ
const q2m = q / m
const μ = 1.0e-20

const E_field = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])
const B_bottle = TP.Field(bottle_B)

"Uniform field: straight field lines, no gradient."
B_uniform(x, t) = SA[0.0, 0.0, 1.0e-8]

"Straight but non-uniform field, so κ = 0 while ∇B ≠ 0."
const GRAD_LENGTH = 2.0
B_gradient(x, t) = SA[0.0, 0.0, 1.0e-8 * (1.0 + x[1] / GRAD_LENGTH)]

"Zero field."
B_zero(x, t) = SA[0.0, 0.0, 0.0]

"Field rotating in z with constant magnitude B0 and wavenumber k."
rotating_B(B0, k) = (x, t) -> SA[B0 * cos(k * x[3]), B0 * sin(k * x[3]), 0.0]

"Fraction of the recorded diagnostics that spent in the full-orbit mode."
fo_frac(sol) = count(==(:FO), sol.stats.adiabaticity.mode) /
    length(sol.stats.adiabaticity.mode)

@testset "Adiabaticity components" begin
    @testset "uniform field" begin
        comps = TP.adiabaticity_components(SA[0.0, 0.0, 0.0], B_uniform, q, m, μ)
        @test comps[1] == 0.0  # ε_curv
        @test comps[2] == 0.0  # ε_gradB
        @test comps[3] == 0.0  # ε_jac (uniform field → zero Jacobian)
    end

    # B = B0 (1 + x1/L) z-hat  ->  κ = 0, ∇B ≠ 0.
    @testset "gradient only" begin
        x = SA[0.5, 0.0, 0.0]
        comps = TP.adiabaticity_components(x, B_gradient, q, m, μ)
        @test comps[1] == 0.0              # straight field lines -> no curvature
        @test comps[2] > 0.0               # finite grad-B length scale
        @test comps[3] ≈ comps[2]          # ε_jac == ε_gradB (κ = 0 → ‖JB‖_F = |∇B|)
    end

    # Curved field: both components generally nonzero.
    @testset "curved field" begin
        x = SA[0.0, 0.0, 1.0]
        comps = TP.adiabaticity_components(x, curved_B_t, q, m, μ)
        @test comps[2] > 0.0               # grad-B present
        @test comps[3] >= max(comps[1], comps[2])  # ε_jac ≥ max(ε_curv, ε_gradB)
    end

    # :curvature selection equals the legacy get_adiabaticity.
    @testset "curvature matches legacy" begin
        for x in (SA[0.0, 0.0, 1.0], SA[1.0, 0.5, 0.0], SA[-0.5, 0.0, 2.0])
            ε_legacy = TP.get_adiabaticity(x, curved_B_t, q, m, μ)
            comps = TP.adiabaticity_components(x, curved_B_t, q, m, μ)
            @test comps[1] ≈ ε_legacy
        end
    end

    # Zero field -> Inf components.
    @testset "zero field" begin
        comps = TP.adiabaticity_components(SA[0.0, 0.0, 0.0], B_zero, q, m, μ)
        @test all(isinf, comps)
    end

    # species convenience method matches explicit q, m.
    @testset "species dispatch" begin
        x = SA[0.0, 0.0, 1.0]
        sp = TP.adiabaticity_components(x, curved_B_t, μ; species = TP.Proton)
        ex = TP.adiabaticity_components(x, curved_B_t, TP.Proton.q, TP.Proton.m, μ)
        @test sp == ex
    end

    # CHIMP jacobian criterion ε_jac for a field that rotates in space with
    # constant magnitude: field lines are straight (κ = 0) and |B| is uniform
    # (∇B = 0), so ε_curv = ε_gradB = 0, yet the B-Jacobian Frobenius norm is
    # finite, giving ε_jac = ρ · k > 0 — exactly the shear/torsion case only
    # CHIMP's criterion captures.
    @testset "rotating field" begin
        B0 = 1.0e-4
        k = 0.5
        x = SA[0.0, 0.0, 0.3]
        comps = TP.adiabaticity_components(x, rotating_B(B0, k), q, m, μ)
        @test comps[1] ≈ 0.0 atol = 1.0e-9  # ε_curv: straight field lines
        @test comps[2] ≈ 0.0 atol = 1.0e-9  # ε_gradB: uniform |B|
        ρ = sqrt(2 * μ * m / B0) / abs(q)
        @test comps[3] ≈ ρ * k rtol = 1.0e-6  # ε_jac = ρ ‖JB‖_F / |B| = ρ k
        @test comps[3] > 0.0
    end
end

@testset "AdaptiveHybrid adiabaticity selection" begin
    x0 = SA[0.0, 0.0, 0.0]
    v0 = SA[5.0e4, 0.0, 1.0e5]
    u0 = vcat(x0, v0)
    T_gyro = 2π / abs(q2m * bottle_B0)
    tspan = (0.0, 30 * T_gyro)
    p = (q2m, m, E_field, B_bottle, TP.ZeroField())

    alg_args = (
        threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro, check_interval = 100,
    )
    prob = TraceHybridProblem(u0, tspan, p)

    # Default (no `adiabaticity`) must equal `:curvature` exactly, i.e. the
    # legacy V1 behaviour is preserved.
    @testset "default is curvature" begin
        sol_def = solve(prob, TP.AdaptiveHybrid(; alg_args...); seed = 1234).u[1]
        sol_curv = solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :curvature); seed = 1234
        ).u[1]
        @test sol_def.retcode == TP.ReturnCode.Success
        @test sol_curv.retcode == TP.ReturnCode.Success
        @test sol_def.t == sol_curv.t
        @test sol_def.u == sol_curv.u
    end

    # All three criteria run and stay accurate/deterministic. With `:both` using
    # OR logic, its full-orbit interval is the union of the curvature and grad-B
    # intervals, so it occupies at least as much full-orbit time as either alone.
    @testset "criteria ordering" begin
        sol_curv = solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :curvature); seed = 1234
        ).u[1]
        sol_grad = solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :gradB); seed = 1234
        ).u[1]
        sol_both = solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :both); seed = 1234
        ).u[1]
        sol_jac = solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :jacobian); seed = 1234
        ).u[1]
        @test sol_curv.retcode == TP.ReturnCode.Success
        @test sol_grad.retcode == TP.ReturnCode.Success
        @test sol_both.retcode == TP.ReturnCode.Success
        @test sol_jac.retcode == TP.ReturnCode.Success

        f_curv = fo_frac(sol_curv)
        f_grad = fo_frac(sol_grad)
        f_both = fo_frac(sol_both)
        f_jac = fo_frac(sol_jac)
        # `:both` is the OR of the other two criteria, so its full-orbit interval
        # is the union and is therefore at least as large as either alone.
        @test f_both >= f_curv
        @test f_both >= f_grad
        # CHIMP's ε_jac dominates ε_curv and ε_gradB, so the jacobian criterion
        # never occupies less full-orbit time than either single criterion.
        @test f_jac >= f_curv
        @test f_jac >= f_grad

        # :curvature reproduces the classic V1 behaviour (mostly guiding-center).
        @test 0.1 < f_curv < 0.5

        # Diagnostics populated for every recorded mode; only the adiabaticity
        # value used by the selected mode is stored (one scalar per check), and
        # it is always non-negative.
        for sol in (sol_curv, sol_grad, sol_both, sol_jac)
            ad = sol.stats.adiabaticity
            @test length(ad.t) == length(ad.components) == length(ad.mode)
            @test length(ad.components) > 0
            @test all(c -> c >= 0.0, ad.components)
        end
        # The initial check point (t = tspan[1]) is recorded identically for all
        # modes (same initial condition, seed = 1234). Verify the per-mode
        # mapping: `:both` is the max of `:curvature`/`:gradB`, and `:jacobian`
        # dominates both single criteria.
        @test sol_both.stats.adiabaticity.components[1] ==
            max(
            sol_curv.stats.adiabaticity.components[1],
            sol_grad.stats.adiabaticity.components[1]
        )
        @test sol_jac.stats.adiabaticity.components[1] >=
            sol_curv.stats.adiabaticity.components[1]
        @test sol_jac.stats.adiabaticity.components[1] >=
            sol_grad.stats.adiabaticity.components[1]
    end

    # A field that rotates in space with constant magnitude: κ = 0 and ∇B = 0,
    # so :curvature and :gradB stay in guiding center, but CHIMP's :jacobian
    # criterion is nonzero and switches to the full orbit.
    @testset "rotating field switches" begin
        B0 = 1.0e-6
        k = 10.0
        B_rot_f = TP.Field(rotating_B(B0, k))
        vperp = sqrt(2 * μ * B0 / m)
        x0 = SA[0.0, 0.0, 0.0]
        u0 = vcat(x0, SA[0.0, 0.0, vperp])
        T_gyro = 2π / abs(q / m * B0)
        tspan = (0.0, 30 * T_gyro)
        p = (q / m, m, E_field, B_rot_f, TP.ZeroField())
        prob_rot = TraceHybridProblem(u0, tspan, p)
        args = (
            threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
            check_interval = 20,
        )

        sol_curv = solve(
            prob_rot, TP.AdaptiveHybrid(; args..., adiabaticity = :curvature);
            seed = 1234
        ).u[1]
        sol_grad = solve(
            prob_rot, TP.AdaptiveHybrid(; args..., adiabaticity = :gradB);
            seed = 1234
        ).u[1]
        sol_jac = solve(
            prob_rot, TP.AdaptiveHybrid(; args..., adiabaticity = :jacobian);
            seed = 1234
        ).u[1]
        @test sol_curv.retcode == TP.ReturnCode.Success
        @test sol_grad.retcode == TP.ReturnCode.Success
        @test sol_jac.retcode == TP.ReturnCode.Success

        f_curv = fo_frac(sol_curv)
        f_grad = fo_frac(sol_grad)
        f_jac = fo_frac(sol_jac)
        # curvature & grad-B see no non-adiabaticity, so they stay in GC.
        @test f_curv < 0.05
        @test f_grad < 0.05
        # the CHIMP jacobian criterion detects the rotating field and switches.
        @test f_jac > max(f_curv, f_grad)
    end

    # Invalid adiabaticity symbol is rejected.
    @testset "invalid criterion" begin
        @test_throws ArgumentError TP.AdaptiveHybrid(; alg_args..., adiabaticity = :bogus)
    end

    # `save_adiabaticity = false` skips the diagnostic buffers entirely while
    # leaving the trajectory unchanged, so users can opt out of the overhead.
    @testset "diagnostics opt out" begin
        sol_save = solve(
            prob, TP.AdaptiveHybrid(; alg_args..., save_adiabaticity = true); seed = 1234
        ).u[1]
        sol_no = solve(
            prob, TP.AdaptiveHybrid(; alg_args..., save_adiabaticity = false); seed = 1234
        ).u[1]
        @test sol_no.retcode == TP.ReturnCode.Success
        @test sol_no.t == sol_save.t
        @test sol_no.u == sol_save.u
        @test sol_no.stats === nothing
    end
end

end # module test_adiabaticity
