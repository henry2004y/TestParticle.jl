using TestParticle
using Test
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using Statistics

@testset "Native GC Solver" begin
    # Dipole field test
    Ek = 5.0e7 # [eV]
    m = TestParticle.mᵢ
    q = TestParticle.qᵢ
    c = TestParticle.c
    Rₑ = TestParticle.Rₑ

    # initial velocity, [m/s]
    v₀ = TestParticle.sph2cart(energy2velocity(Ek; q, m), π / 4, 0.0)
    # initial position, [m]
    r₀ = TestParticle.sph2cart(2.5 * Rₑ, π / 2, 0.0)
    stateinit = [r₀..., v₀...]
    tspan = (0.0, 1.0)

    # 1. Prepare GC parameters
    stateinit_gc, param_gc = TestParticle.prepare_gc(
        stateinit, TestParticle.ZeroField(), TestParticle.getB_dipole;
        species = Proton
    )

    # 2. Solve using DiffEq
    prob_diffeq = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
    sol_diffeq = solve(prob_diffeq, Vern9(); reltol = 1.0e-8, abstol = 1.0e-8)

    # 3. Solve using Native RK4
    prob_native = TraceGCProblem(stateinit_gc, tspan, param_gc)

    dt = 1.0e-4 # Reasonable step size for RK4
    sol_native = TestParticle.solve(prob_native; dt = dt, alg = :rk4)

    # Compare end point
    u_diffeq = sol_diffeq[end]
    u_native = sol_native[1].u[end]

    # Position difference
    @test norm(u_diffeq[1:3] - u_native[1:3]) / norm(u_diffeq[1:3]) < 1.0e-3

    # Parallel velocity difference
    @test abs(u_diffeq[4] - u_native[4]) / abs(u_diffeq[4]) < 1.0e-3

end
