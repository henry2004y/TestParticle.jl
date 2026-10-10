if !isdefined(Main, :test_common)
    include("test_common.jl")
end

module test_time_dependent_fields

using Test
using TestParticle
import TestParticle as TP
using OrdinaryDiffEq
using StaticArrays
using ..test_common: uniform_grid, cartesian_grid

"Oscillating electric field along +x."
E_osc(r, t) = SA[5.0e-11 * sin(2π * t), 0, 0]  # [V/m]

"Oscillating magnetic field along +z."
B_osc(r, t) = SA[0, 0, 1.0e-11 * cos(2π * t)]  # [T]

"Constant force field along +y."
F_const(r) = SA[0, 9.10938356e-42, 0]  # [N]

@testset "time-dependent fields" begin
    x, y, z = uniform_grid(15, 20, 25)
    F = fill(0.0, 3, length(x), length(y), length(z)) # [N]
    F[2, :, :, :] .= 9.10938356e-42

    grid = cartesian_grid(x, y, z)

    x0 = [0.0, 0.0, 0.0] # initial position, [m]
    u0 = [0.0, 1.0, 1.0] # initial velocity, [m/s]
    stateinit = [x0..., u0...]
    tspan = (0.0, 1.0)

    param = prepare(grid, E_osc, B_osc, F; species = Electron)
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1, 2, 3])
    @test isapprox(sol[1, end], -1.2828663442681638, rtol = 1.0e-2) &&
        isapprox(sol[2, end], 1.5780464321537067, rtol = 1.0e-2) &&
        isapprox(sol[3, end], 1.0, rtol = 1.0e-2)
    @test get_gc([stateinit..., 0.0], param)[1] ≈ -0.5685630103565723

    param = prepare(E_osc, B_osc, F_const; species = Electron)
    _, _, _, _, F_out = param
    @test F_out(x0)[2] ≈ 9.10938356e-42

    stateinit = SA[x0..., u0...]

    prob = ODEProblem(trace, stateinit, tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1, 2, 3])
    @test isapprox(sol[1, end], -1.2828663442681638, rtol = 1.0e-2) &&
        isapprox(sol[2, end], 1.5780464321537067, rtol = 1.0e-2) &&
        isapprox(sol[3, end], 1.0, rtol = 1.0e-2)
end

end # module test_time_dependent_fields
