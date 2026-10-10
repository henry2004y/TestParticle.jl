if !isdefined(Main, :test_common)
    include("test_common.jl")
end

module test_relativistic

using Test
using TestParticle
import TestParticle as TP
using OrdinaryDiffEq
using StaticArrays
using ..test_common: uniform_B_strong, zero_E

"Uniform unit magnetic field along +z, for the dimensionless form."
B_unit(x) = SA[0.0, 0.0, 1.0]

"Potential well over the square 0 ≤ x, y ≤ 100 m, holding 1.0e5 V/m."
function E_well(xu)
    Ex = 0 <= xu[1] <= 100 ? -1.0e5 : 0.0
    Ey = 0 <= xu[2] <= 100 ? -1.0e5 : 0.0
    return SA[Ex, Ey, 0.0]
end

@testset "relativistic particle" begin
    x0 = [10.0, 10.0, 0.0] # initial position, [m]
    u0 = [0.0, 0.0, 0.0] # initial velocity γv, [m/s]
    tspan = (0.0, 2.0e-7)
    stateinit = [x0..., u0...]

    param = prepare(E_well, uniform_B_strong; species = Electron)
    prob = ODEProblem(trace_relativistic!, stateinit, tspan, param)
    sol = solve(prob, Vern6(); dt = 1.0e-10, adaptive = false)

    # The kinetic energy [eV] of the electron has to match the electric
    # potential energy gained along the way.
    x = sol.u[end][1:3]
    Δx = x[1] - x0[1] + x[2] - x0[2]
    KE_end = get_energy(sol[4:6, end]; m = TP.mₑ, q = TP.qₑ)
    @test isapprox(KE_end / Δx, 1.0e5, rtol = 1.0e-2)
    KE_all = get_energy(sol; m = TP.mₑ, q = TP.qₑ)
    @test isapprox(KE_all[end] / Δx, 1.0e5, rtol = 1.0e-2)
    # Convert to velocity
    @test isapprox(get_velocity(sol)[1, end], -1.6588724708347203e7, rtol = 1.0e-2)

    prob = ODEProblem(trace_relativistic, SA[stateinit...], tspan, param)
    sol = solve(prob, Vern6(); dt = 1.0e-10, adaptive = false)

    x = sol[1:3, end]
    Δx = x[1] - x0[1] + x[2] - x0[2]
    KE_end = get_energy(sol[4:6, end]; m = TP.mₑ, q = TP.qₑ)
    @test isapprox(KE_end / Δx, 1.0e5, rtol = 1.0e-2)

    # Tracing relativistic particles in dimensionless units
    param = prepare(zero_E, B_unit; m = 1, q = 1)
    tspan = (0.0, 1.0) # 1/2π period
    stateinit = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0]
    prob = ODEProblem(trace_relativistic_normalized!, stateinit, tspan, param)
    sol = solve(prob, Vern6())
    @test isapprox(sol.u[end][1], 0.38992532495827226, rtol = 1.0e-2)
    stateinit = zeros(6)
    prob = ODEProblem(trace_relativistic_normalized!, stateinit, tspan, param)
    sol = solve(prob, Vern6(); abstol = 1.0e-8, reltol = 1.0e-8, dt = 1.0e-4)
    @test isapprox(sol[1, end], 0.0, atol = 1.0e-2) && length(sol.u) == 3

    stateinit = SA[0.0, 0.0, 0.0, 0.5, 0.0, 0.0]
    prob = ODEProblem(trace_relativistic_normalized, stateinit, tspan, param)
    sol = solve(prob, Vern6())
    @test isapprox(sol[1, end], 0.38992532495827226, rtol = 1.0e-2)

    stateinit = @SVector zeros(6)
    prob = ODEProblem(trace_relativistic_normalized, stateinit, tspan, param)
    sol = solve(prob, Vern6(); abstol = 1.0e-8, reltol = 1.0e-8, dt = 1.0e-4)
    @test isapprox(sol[1, end], 0.0, atol = 1.0e-2) && length(sol.u) == 3
end

end # module test_relativistic
