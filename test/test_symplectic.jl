module test_symplectic

using Test
using TestParticle
using OrdinaryDiffEq
using OrdinaryDiffEqSymplecticRK
using OrdinaryDiffEqSDIRK
using StaticArrays
using LinearAlgebra: norm

"Kinetic energy of the velocity half of a `DynamicalODEProblem` state."
kinetic_energy(v; m) = 0.5 * m * norm(v)^2

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

        # H = K + V with V = -q E_y y² / 2
        get_H(u) = kinetic_energy(u.x[1]; m) - 0.5 * q * E_y * u.x[2][2]^2

        @test get_H(sol.u[1]) ≈ get_H(sol.u[end]) rtol = 1.0e-4
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

        # McAte2 is an explicit symplectic integrator designed for separable
        # Hamiltonians (H = T(p) + V(q)). The Lorentz force Hamiltonian
        # H = |p-qA|² / 2m is non-separable when B ≠ 0 due to the vector
        # potential A(x). Thus, explicit partitioned Runge-Kutta methods like
        # McAte2 are not exactly symplectic for the Lorentz force and may
        # exhibit energy drift, requiring smaller dt or relaxed tolerance.
        sol = solve(prob, McAte2(), dt = 0.01, adaptive = false)

        K0 = kinetic_energy(v0; m)
        @test K0 ≈ kinetic_energy(sol.u[end].x[1]; m) rtol = 0.04

        # ImplicitMidpoint is symplectic for all Hamiltonian systems, including
        # non-separable ones.
        prob_ode = ODEProblem(trace_normalized!, [x0..., v0...], tspan, param)
        sol_im = solve(prob_ode, ImplicitMidpoint(), dt = 0.01, adaptive = false)
        @test K0 ≈ kinetic_energy(sol_im.u[end][4:6]; m) rtol = 1.0e-12
    end
end

end # module test_symplectic
