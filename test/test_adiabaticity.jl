using TestParticle
import TestParticle as TP
using StaticArrays
using LinearAlgebra
using Test

@testset "Adiabaticity components" begin
    m = TP.mᵢ
    q = TP.qᵢ
    μ = 1.0e-20

    # Uniform field: straight field lines, no gradient -> all zero.
    let
        B(x, t) = SA[0.0, 0.0, 1.0e-8]
        comps = TP.adiabaticity_components(SA[0.0, 0.0, 0.0], B, q, m, μ)
        @test comps[1] == 0.0  # ε_curv
        @test comps[2] == 0.0  # ε_gradB
        @test comps[3] == 0.0  # ε_jac (uniform field → zero Jacobian)
    end

    # Straight but non-uniform field (grad-B only, no curvature):
    # B = B0 (1 + x1/L) z-hat  ->  κ = 0, ∇B ≠ 0.
    let
        L = 2.0
        B(x, t) = SA[0.0, 0.0, 1.0e-8 * (1.0 + x[1] / L)]
        x = SA[0.5, 0.0, 0.0]
        comps = TP.adiabaticity_components(x, B, q, m, μ)
        @test comps[1] == 0.0              # straight field lines -> no curvature
        @test comps[2] > 0.0               # finite grad-B length scale
        @test comps[3] ≈ comps[2]          # ε_jac == ε_gradB (κ = 0 → ‖JB‖_F = |∇B|)
    end

    # Curved field: both components generally nonzero.
    let
        function B(x, t)
            θ = atan(x[3] / (x[1] + 3))
            r = sqrt((x[1] + 3)^2 + x[3]^2)
            return SA[-1.0e-6 * sin(θ) / r, 0.0, 1.0e-6 * cos(θ) / r]
        end
        x = SA[0.0, 0.0, 1.0]
        comps = TP.adiabaticity_components(x, B, q, m, μ)
        @test comps[2] > 0.0               # grad-B present
        @test comps[3] >= max(comps[1], comps[2])  # ε_jac ≥ max(ε_curv, ε_gradB)
    end

    # :curvature selection equals the legacy get_adiabaticity.
    let
        function B(x, t)
            θ = atan(x[3] / (x[1] + 3))
            r = sqrt((x[1] + 3)^2 + x[3]^2)
            return SA[-1.0e-6 * sin(θ) / r, 0.0, 1.0e-6 * cos(θ) / r]
        end
        for x in (SA[0.0, 0.0, 1.0], SA[1.0, 0.5, 0.0], SA[-0.5, 0.0, 2.0])
            ε_legacy = TP.get_adiabaticity(x, B, q, m, μ)
            comps = TP.adiabaticity_components(x, B, q, m, μ)
            @test comps[1] ≈ ε_legacy
        end
    end

    # Zero field -> Inf components.
    let
        B(x, t) = SA[0.0, 0.0, 0.0]
        comps = TP.adiabaticity_components(SA[0.0, 0.0, 0.0], B, q, m, μ)
        @test all(isinf, comps)
    end

    # species convenience method matches explicit q, m.
    let
        function B(x, t)
            θ = atan(x[3] / (x[1] + 3))
            r = sqrt((x[1] + 3)^2 + x[3]^2)
            return SA[-1.0e-6 * sin(θ) / r, 0.0, 1.0e-6 * cos(θ) / r]
        end
        x = SA[0.0, 0.0, 1.0]
        sp = TP.adiabaticity_components(x, B, μ; species = TP.Proton)
        ex = TP.adiabaticity_components(x, B, TP.Proton.q, TP.Proton.m, μ)
        @test sp == ex
    end

    # CHIMP jacobian criterion ε_jac for a field that rotates in space with
    # constant magnitude: field lines are straight (κ = 0) and |B| is uniform
    # (∇B = 0), so ε_curv = ε_gradB = 0, yet the B-Jacobian Frobenius norm is
    # finite, giving ε_jac = ρ · k > 0 — exactly the shear/torsion case only
    # CHIMP's criterion captures.
    let
        B0 = 1.0e-4
        k = 0.5
        B(x, t) = SA[B0 * cos(k * x[3]), B0 * sin(k * x[3]), 0.0]
        x = SA[0.0, 0.0, 0.3]
        comps = TP.adiabaticity_components(x, B, q, m, μ)
        @test comps[1] ≈ 0.0 atol = 1.0e-9  # ε_curv: straight field lines
        @test comps[2] ≈ 0.0 atol = 1.0e-9  # ε_gradB: uniform |B|
        ρ = sqrt(2 * μ * m / B0) / abs(q)
        @test comps[3] ≈ ρ * k rtol = 1.0e-6  # ε_jac = ρ ‖JB‖_F / |B| = ρ k
        @test comps[3] > 0.0
    end
end

@testset "AdaptiveHybrid adiabaticity selection" begin
    m = TP.mᵢ
    q = TP.qᵢ
    q2m = q / m
    E_field = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

    # Bottle field with both curvature and grad-B.
    B0 = 1.0e-4
    α = 1.0e-2
    function bottle_B(x, t)
        Bz = B0 * (1 + α * x[3]^2)
        Bx = -B0 * α * x[1] * x[3]
        By = -B0 * α * x[2] * x[3]
        return SA[Bx, By, Bz]
    end
    B_bottle = TP.Field(bottle_B)

    x0 = SA[0.0, 0.0, 0.0]
    v0 = SA[5.0e4, 0.0, 1.0e5]
    u0 = vcat(x0, v0)
    Ω = abs(q2m) * B0
    T_gyro = 2π / Ω
    tspan = (0.0, 30 * T_gyro)
    p = (q2m, m, E_field, B_bottle, TP.ZeroField())

    alg_args = (
        threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro, check_interval = 100
    )
    prob = TraceHybridProblem(u0, tspan, p)

    # Helper: FO-mode fraction recorded in the diagnostics.
    fo_frac(sol) = count(==(:FO), sol.stats.adiabaticity.mode) /
                   length(sol.stats.adiabaticity.mode)

    # Default (no `adiabaticity`) must equal `:curvature` exactly, i.e. the
    # legacy V1 behaviour is preserved.
    let
        sol_def = TP.solve(prob, TP.AdaptiveHybrid(; alg_args...); seed = 1234).u[1]
        sol_curv = TP.solve(
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
    let
        sol_curv = TP.solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :curvature); seed = 1234
        ).u[1]
        sol_grad = TP.solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :gradB); seed = 1234
        ).u[1]
        sol_both = TP.solve(
            prob, TP.AdaptiveHybrid(; alg_args..., adiabaticity = :both); seed = 1234
        ).u[1]
        sol_jac = TP.solve(
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
              max(sol_curv.stats.adiabaticity.components[1],
                  sol_grad.stats.adiabaticity.components[1])
        @test sol_jac.stats.adiabaticity.components[1] >=
              sol_curv.stats.adiabaticity.components[1]
        @test sol_jac.stats.adiabaticity.components[1] >=
              sol_grad.stats.adiabaticity.components[1]
    end

    # A field that rotates in space with constant magnitude: κ = 0 and ∇B = 0,
    # so :curvature and :gradB stay in guiding center, but CHIMP's :jacobian
    # criterion is nonzero and switches to the full orbit.
    let
        μ = 1.0e-20
        B0 = 1.0e-6
        k = 10.0
        B_rot(x, t) = SA[B0 * cos(k * x[3]), B0 * sin(k * x[3]), 0.0]
        B_rot_f = TP.Field(B_rot)
        vperp = sqrt(2 * μ * B0 / m)
        x0 = SA[0.0, 0.0, 0.0]
        v0 = SA[0.0, 0.0, vperp]
        u0 = vcat(x0, v0)
        Ω = abs(q / m) * B0
        T_gyro = 2π / Ω
        tspan = (0.0, 30 * T_gyro)
        p = (q / m, m, E_field, B_rot_f, TP.ZeroField())
        prob_rot = TraceHybridProblem(u0, tspan, p)
        args = (threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
                check_interval = 20)

        sol_curv = TP.solve(
            prob_rot, TP.AdaptiveHybrid(; args..., adiabaticity = :curvature);
            seed = 1234
        ).u[1]
        sol_grad = TP.solve(
            prob_rot, TP.AdaptiveHybrid(; args..., adiabaticity = :gradB);
            seed = 1234
        ).u[1]
        sol_jac = TP.solve(
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
    @test_throws ArgumentError TP.AdaptiveHybrid(; alg_args..., adiabaticity = :bogus)

    # `save_adiabaticity = false` skips the diagnostic buffers entirely while
    # leaving the trajectory unchanged, so users can opt out of the overhead.
    let
        sol_save = TP.solve(
            prob, TP.AdaptiveHybrid(; alg_args..., save_adiabaticity = true); seed = 1234
        ).u[1]
        sol_no = TP.solve(
            prob, TP.AdaptiveHybrid(; alg_args..., save_adiabaticity = false); seed = 1234
        ).u[1]
        @test sol_no.retcode == TP.ReturnCode.Success
        @test sol_no.t == sol_save.t
        @test sol_no.u == sol_save.u
        @test sol_no.stats === nothing
    end
end
