using Test
using TestParticle
using OrdinaryDiffEq
using OrdinaryDiffEqSymplecticRK
using OrdinaryDiffEqSDIRK
using StaticArrays
using LinearAlgebra: norm

@testset "Symplectic Solvers" begin
    @testset "Separable E-field" begin
        q = 1.0
        m = 1.0
        E_y = 1.0
        v0 = [0.0, 1.0, 0.0]
        x0 = [1.0, 0.0, 0.0]
        tspan = (0.0, 1.0)

        # Spatially linear E field in Y direction: E = E_y * y * y_hat
        # H = p^2/2m - q * E_y * y^2 / 2 (Separable)
        spatial_linear_E(x, t) = SA[0.0, E_y * x[2], 0.0]
        param = prepare(spatial_linear_E, ZeroField(); q, m)

        # get_dv! and get_dx! are exported from TestParticle
        prob = DynamicalODEProblem(get_dv!, get_dx!, v0, x0, tspan, param)
        # McAte2 is an explicit symplectic integrator. For this separable system,
        # it should conserve energy well with O(dt^2) error.
        sol = solve(prob, McAte2(), dt = 0.01, adaptive = false)

        # Energy conservation check (Total Energy H = K + V)
        function get_H(u)
            v = u.x[1]
            x = u.x[2]
            K = 0.5 * m * norm(v)^2
            V = -0.5 * q * E_y * x[2]^2
            return K + V
        end

        H0 = get_H(sol.u[1])
        H_final = get_H(sol.u[end])
        @test H0 ≈ H_final rtol = 1.0e-4
    end

    @testset "Constant B-field" begin
        q = 1.0
        m = 1.0
        B_z = 1.0
        v0 = [1.0, 0.0, 0.0]
        x0 = [0.0, 0.0, 0.0]
        tspan = (0.0, 2π)

        constant_B(x, t) = SA[0.0, 0.0, B_z]
        param = prepare(ZeroField(), constant_B; q, m)

        prob = DynamicalODEProblem(get_dv!, get_dx!, v0, x0, tspan, param)

        # McAte2 is an explicit symplectic integrator designed for separable Hamiltonians (H = T(p) + V(q)).
        # The Lorentz force Hamiltonian H = |p-qA|^2 / 2m is non-separable when B != 0 due to the vector potential A(x).
        # Thus, explicit partitioned Runge-Kutta methods like McAte2 are not exactly symplectic for the Lorentz force
        # and may exhibit energy drift, requiring smaller dt or relaxed tolerance.
        sol = solve(prob, McAte2(), dt = 0.01, adaptive = false)

        K0 = 0.5 * m * norm(v0)^2
        v_final = sol.u[end].x[1]
        K_final = 0.5 * m * norm(v_final)^2
        @test K0 ≈ K_final rtol = 0.04

        # ImplicitMidpoint is symplectic for all Hamiltonian systems, including non-separable ones.
        u0 = [x0..., v0...]
        prob_ode = ODEProblem(trace_normalized!, u0, tspan, param)
        sol_im = solve(prob_ode, ImplicitMidpoint(), dt = 0.01, adaptive = false)
        K_im = 0.5 * m * norm(sol_im.u[end][4:6])^2
        @test K0 ≈ K_im rtol = 1.0e-12
    end
end
