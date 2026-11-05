module test_gc_relativistic

using TestParticle
using TestParticle: c, Rₑ, BMoment_Earth, full_to_gc, gc_to_full, get_gc,
    get_gc_velocity, get_energy, get_adiabaticity, get_work_channels, sph2cart
using StaticArrays
using Test
using LinearAlgebra
using SciMLBase
using OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqVerner
import Magnetostatics as MS

const dipole = MS.Dipole(BMoment_Earth)
dipole_field(r) = dipole(r)

@testset "RelativisticGuidingCenter" begin
    @testset "FO <-> GC conversion" begin
        let B0 = 1.0, E0 = 1.0e5
            B_u(r) = SA[0.0, 0.0, B0]
            E_u(r) = SA[0.0, E0, 0.0]
            param_B = prepare(; B = B_u, species = Electron)
            param_EB = prepare(; E = E_u, B = B_u, species = Electron)

            # Pure B: exact round-trip GC -> FO -> GC at relativistic energy (γ ≈ 3)
            R = SA[1.0, -2.0, 0.5]
            upar = 2.0 * c
            u_perp = 2.0 * c
            μ = Electron.m * u_perp^2 / (2 * B0)
            phase = π / 3
            gc_state = SA[R[1], R[2], R[3], upar, μ]

            xu = gc_to_full(gc_state, param_B; phase, model = RelativisticGuidingCenter())
            @test length(xu) == 6
            # Proper velocity magnitude matches √(u_∥² + u_⊥²)
            @test norm(xu[4:6]) ≈ hypot(upar, u_perp) rtol = 1.0e-12

            @test full_to_gc(xu, param_B; model = RelativisticGuidingCenter()) ≈ gc_state rtol = 1.0e-12

            # With E × B drift (v_E = E0 / B0 = 1e5 m/s along +x), μ = 0
            gc_drift = SA[0.0, 0.0, 0.0, upar, 0.0]
            γ0 = sqrt(1 + (upar / c)^2)
            xu_eb = gc_to_full(gc_drift, param_EB; model = RelativisticGuidingCenter())
            @test xu_eb[4] ≈ γ0 * (E0 / B0) rtol = 1.0e-12
            @test xu_eb[5] ≈ 0.0 atol = 1.0e-12
            @test xu_eb[6] ≈ upar rtol = 1.0e-12

            gc_eb_back = full_to_gc(xu_eb, param_EB; model = RelativisticGuidingCenter())
            @test gc_eb_back[1:3] ≈ SA[0.0, 0.0, 0.0] atol = 1.0e-12
            @test gc_eb_back[4] ≈ upar rtol = 1.0e-12
            @test gc_eb_back[5] ≈ 0.0 atol = 1.0e-30

            # get_gc with use_vE = true subtracts γ v_E for relativistic states
            X_gc = get_gc(
                xu_eb, param_EB; use_vE = true, model = RelativisticGuidingCenter()
            )
            # In uniform B with μ = 0, subtracting γ v_E leaves zero perp velocity -> X = x
            @test X_gc ≈ xu_eb[1:3] atol = 1.0e-6
        end
    end

    @testset "Non-relativistic limit matches GuidingCenter" begin
        let B_grad(r) = SA[0.0, 0.0, 1.0e-5 * (1 + 1.0e-6 * r[1])]
            v0 = SA[1.0e4, 0.0, 2.0e4] # v / c ~ 7e-5
            x0 = SA[0.0, 0.0, 0.0]
            u0 = SA[x0..., v0...]

            p = prepare(; B = B_grad, species = Proton)

            tspan = (0.0, 1.0)
            sol_nr = solve(
                TraceGCProblem(u0, tspan, p), Vern7();
                reltol = 1.0e-10, abstol = 1.0e-10
            )
            sol_rel = solve(
                TraceProblem{RelativisticGuidingCenter}(u0, tspan, p), Vern7();
                reltol = 1.0e-10, abstol = 1.0e-10
            )
            @test sol_rel.retcode == ReturnCode.Success
            @test sol_rel.u[end] ≈ sol_nr.u[end] rtol = 1.0e-6
        end
    end

    @testset "Relativistic electron in dipole field vs RelativisticOrbit" begin
        let
            # 2 MeV electron at L = 4 Rₑ, equatorial pitch angle 60°
            γ0 = 1 + 2.0e6 / 5.1099895e5 # ≈ 4.914
            v_mag = c * sqrt(1 - inv(γ0^2))
            α = π / 3
            u_par = γ0 * v_mag * cos(α)
            u_perp = γ0 * v_mag * sin(α)
            r0 = SA[4.0 * Rₑ, 0.0, 0.0]
            # At equator (x > 0, y = z = 0), dipole B points in -z
            u0_fo = SA[r0[1], r0[2], r0[3], u_perp, 0.0, u_par]

            tspan = (0.0, 0.5)

            p_fo = prepare(; B = dipole_field, species = Electron)
            prob_fo = TraceProblem{RelativisticOrbit}(u0_fo, tspan, p_fo)
            sol_fo = solve(prob_fo, Vern7(); reltol = 1.0e-9, abstol = 1.0e-9)

            p_gc = p_fo
            prob_gc = TraceProblem{RelativisticGuidingCenter}(u0_fo, tspan, p_gc)
            u0_gc = prob_gc.u0
            sol_gc = solve(prob_gc, Vern7(); reltol = 1.0e-9, abstol = 1.0e-9)

            @test sol_gc.retcode == ReturnCode.Success
            @test length(sol_gc.t) < length(sol_fo.t)

            # Guiding-center trajectory tracks full-orbit guiding center
            X_fo_end = get_gc(sol_fo.u[end], p_fo)
            X_gc_end = sol_gc.u[end][1:3]
            @test norm(X_gc_end - X_fo_end) / (4.0 * Rₑ) < 1.0e-3

            # Energy conservation in static dipole field
            E_gc = get_energy(sol_gc)
            E_fo = get_energy(sol_fo)
            @test E_gc[1] ≈ E_fo[1] rtol = 1.0e-3
            @test maximum(abs.(E_gc .- E_gc[1])) / E_gc[1] < 1.0e-7

            # In-place and out-of-place equations agree
            du = zeros(5)
            trace_gc_relativistic!(du, u0_gc, p_gc, 0.0)
            @test du ≈ trace_gc_relativistic(u0_gc, p_gc, 0.0) rtol = 1.0e-14
            @test get_gc_velocity(
                u0_gc, p_gc, 0.0; model = RelativisticGuidingCenter()
            ) ≈ du[1:3] rtol = 1.0e-14

            # Diagnostics: adiabaticity and work rates
            ϵ_gc = get_adiabaticity(sol_gc)
            ϵ_fo = get_adiabaticity(sol_fo)
            @test ϵ_gc[1] ≈ ϵ_fo[1] rtol = 1.0e-2
            @test all(isfinite, ϵ_gc)

            w_gc = get_work_rates(
                sol_gc.u[1], p_gc, 0.0; model = RelativisticGuidingCenter()
            )
            @test all(iszero, w_gc)
            ch_gc = get_work_channels(
                sol_gc.u[1], p_gc, 0.0; model = RelativisticGuidingCenter()
            )
            @test all(iszero, ch_gc)
        end
    end

    @testset "Normalized units (q = 1, m = 1, c = 1)" begin
        let
            # Magnetic mirror B_z(z) = 1 + 0.1 z^2, B_r = -0.1 r z (∇⋅B = 0)
            B_mirror(r) = SA[-0.1 * r[1] * r[3], -0.1 * r[2] * r[3], 1.0 + 0.1 * r[3]^2]
            E_z(r) = SA[0.0, 0.0, 0.05]

            # Pure mirror: energy conserved, parallel momentum mirrors
            u0_fo = SA[0.0, 0.0, 0.0, 1.5, 0.0, 1.0] # γ = √(1 + 3.25) ≈ 2.06
            p_gc = prepare(; B = B_mirror, q = 1.0, m = 1.0, c = 1.0)
            prob = TraceProblem{RelativisticGuidingCenter}(u0_fo, (0.0, 10.0), p_gc)
            sol = solve(prob, Vern7(); reltol = 1.0e-10, abstol = 1.0e-10)
            @test sol.retcode == ReturnCode.Success
            E = get_energy(sol)
            @test maximum(abs.(E .- E[1])) / E[1] < 1.0e-8

            # With uniform E_z along uniform B_z: du_∥/dt = (q/m) E_z, dH/dt = 0
            B_z(r) = SA[0.0, 0.0, 1.0]
            p_e = prepare(; E = E_z, B = B_z, q = 1.0, m = 1.0, c = 1.0)
            prob_e = TraceProblem{RelativisticGuidingCenter}(
                SA[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], (0.0, 2.0), p_e
            )
            sol_e = solve(prob_e, RK4(); dt = 1.0e-2, adaptive = false)
            @test sol_e.retcode == ReturnCode.Success
            @test sol_e.u[end][4] ≈ 0.05 * 2.0 rtol = 1.0e-8

            # Work rate P_par = q (u_∥ / γ) E_z
            u_end = sol_e.u[end]
            γ_end = sqrt(1 + u_end[4]^2 + 2 * u_end[5] * 1.0)
            rates = get_work_rates(u_end, p_e, 2.0; model = RelativisticGuidingCenter())
            @test rates[1] ≈ (u_end[4] / γ_end) * 0.05 rtol = 1.0e-12
        end
    end

    @testset "Problem construction and Float32 preservation" begin
        let B_u(r) = SA[0.0f0, 0.0f0, 1.0f-5]
            u0 = SA[0.0f0, 0.0f0, 0.0f0, 1.0f8, 0.0f0, 1.0f8]
            p_gc = prepare(; B = B_u, species = Electron)
            prob = TraceProblem{RelativisticGuidingCenter}(u0, (0.0f0, 1.0f-4), p_gc)
            u0_gc = prob.u0
            @test u0_gc isa SVector{5, Float32}
            sol = solve(prob, Vern7())
            @test sol.retcode == ReturnCode.Success
            @test eltype(sol.u[end]) === Float32

            # Constructing from equations directly infers RelativisticGuidingCenter
            prob_from_f = TraceProblem(
                trace_gc_relativistic, u0_gc, (0.0, 1.0e-4), p_gc
            )
            @test prob_from_f isa TraceProblem{RelativisticGuidingCenter}

            # Mismatched equations throw ArgumentError
            @test_throws ArgumentError TraceProblem{RelativisticGuidingCenter}(
                trace_gc!, u0_gc, (0.0, 1.0e-4), p_gc
            )
        end
    end
end

end
