module test_canonical

using Test
using TestParticle
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using StaticArrays
using LinearAlgebra: norm

@testset "Canonical Form" begin
    @testset "PotentialField Construction & AD" begin
        let
            # Analytical scalar and vector potentials
            phi(x, t) = x[1]^2 + 2 * x[2]^2
            A(x, t) = SA[0.0, 3.0 * x[1], 0.0]

            pf = PotentialField(phi, A)

            x0 = SA[1.0, 2.0, 0.0]
            # ∇ϕ = [2x₁, 4x₂, 0] = [2, 8, 0]
            @test pf.grad_phi(x0, 0.0) ≈ SA[2.0, 8.0, 0.0]
            # ∇A: (∂A_j / ∂x_i) -> ∂A_y / ∂x = 3.0 at row 1, col 2
            expected_grad_A = SA[0.0 3.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
            @test pf.grad_A(x0, 0.0) ≈ expected_grad_A

            # PotentialField with ZeroField
            pf_zero_phi = PotentialField(A)
            @test pf_zero_phi.grad_phi isa ZeroField
            @test pf_zero_phi.grad_A(x0, 0.0) ≈ expected_grad_A

            pf_zero_A = PotentialField(phi, ZeroField())
            @test pf_zero_A.grad_A(x0, 0.0) == zero(SMatrix{3, 3, Float64, 9})
        end
    end

    @testset "Non-Relativistic Landau Gauge" begin
        let
            q = 1.0
            m = 1.0
            B0 = 1.0
            Ω = q * B0 / m
            T_period = 2π / Ω

            A_landau(x) = SA[0.0, B0 * x[1], 0.0]
            pf = PotentialField(A_landau)
            param = prepare(pf; q = q, m = m)

            x0 = SA[0.0, 0.0, 0.0]
            v0 = SA[1.0, 0.0, 0.0]
            u0 = velocity_to_canonical(x0, v0, param; relativistic = false)

            # In Landau gauge, p_y = m vy + q A_y = 0 initially
            @test u0[5] ≈ 0.0 atol = 1.0e-15

            tspan = (0.0, 5 * T_period)
            dt = T_period / 40

            prob = TraceCanonicalProblem(u0, tspan, param; relativistic = false)
            sol = solve(prob, ImplicitMidpoint(); dt = dt, adaptive = false)

            # 1. Exact conservation of cyclic momentum py to machine precision
            py_start = sol.u[1][5]
            py_end = sol.u[end][5]
            @test py_start ≈ py_end atol = 1.0e-14

            # 2. Energy conservation to machine precision with ImplicitMidpoint
            H_start = canonical_hamiltonian(sol.u[1], param; relativistic = false)
            H_end = canonical_hamiltonian(sol.u[end], param; relativistic = false)
            @test H_start ≈ H_end rtol = 1.0e-12

            # 3. Trajectory and velocity conversion
            v_end = canonical_to_velocity(sol.u[end], param; relativistic = false)
            @test norm(v_end) ≈ norm(v0) rtol = 1.0e-12

            # In-place mutable vector support
            u0_mut = [x0..., (u0[4:6])...]
            prob_mut = TraceCanonicalProblem(u0_mut, tspan, param; relativistic = false)
            sol_mut = solve(prob_mut, ImplicitMidpoint(); dt = dt, adaptive = false)
            @test sol_mut.u[end][5] ≈ py_start atol = 1.0e-14
        end
    end

    @testset "Relativistic Landau Gauge" begin
        let
            q = 1.0
            m = 1.0
            c_val = 1.0
            B0 = 1.0
            Ω0 = q * B0 / m
            T_period = 2π / Ω0

            A_landau(x) = SA[0.0, B0 * x[1], 0.0]
            pf = PotentialField(A_landau)
            param = prepare(pf; q = q, m = m, c = c_val)

            x0 = SA[0.0, 0.0, 0.0]
            v0 = SA[0.6 * c_val, 0.0, 0.0]
            u0 = velocity_to_canonical(x0, v0, param; relativistic = true)

            tspan = (0.0, 5 * T_period)
            dt = T_period / 50

            prob = TraceCanonicalProblem(u0, tspan, param; relativistic = true)
            sol = solve(prob, ImplicitMidpoint(); dt = dt, adaptive = false)

            # Cyclic momentum py must be conserved to machine precision
            @test sol.u[1][5] ≈ sol.u[end][5] atol = 1.0e-14

            # Symplectic energy error is bounded
            H_start = canonical_hamiltonian(sol.u[1], param; relativistic = true)
            H_end = canonical_hamiltonian(sol.u[end], param; relativistic = true)
            @test (abs(H_end - H_start) / H_start) < 1.0e-3

            # Velocity recovery
            v_end = canonical_to_velocity(sol.u[end], param; relativistic = true)
            @test norm(v_end) ≈ norm(v0) rtol = 1.0e-3
        end
    end

    @testset "Electrostatic Harmonic Oscillator" begin
        let
            q = 1.0
            m = 1.0
            k = 4.0
            phi_harmonic(x, t) = 0.5 * k * x[1]^2

            param = prepare_canonical(phi_harmonic, ZeroField(); q = q, m = m)

            x0 = SA[1.0, 0.0, 0.0]
            v0 = SA[0.0, 0.0, 0.0]
            u0 = velocity_to_canonical(x0, v0, param; relativistic = false)

            ω = √(k * q / m) # ω = 2.0
            T_osc = 2π / ω
            tspan = (0.0, 2 * T_osc)
            dt = T_osc / 100

            prob = TraceCanonicalProblem(u0, tspan, param; relativistic = false)
            sol = solve(prob, ImplicitMidpoint(); dt = dt, adaptive = false)

            # Energy should be conserved
            H_start = canonical_hamiltonian(sol.u[1], param; relativistic = false)
            H_end = canonical_hamiltonian(sol.u[end], param; relativistic = false)
            @test H_start ≈ H_end rtol = 1.0e-12

            # After integer periods, particle returns to x0
            @test sol.u[end][1:3] ≈ x0 rtol = 1.0e-3
        end
    end
end

end # module test_canonical
