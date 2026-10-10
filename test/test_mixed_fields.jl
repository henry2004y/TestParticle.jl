module test_mixed_fields

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using Test
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using ..test_common: uniform_grid, cartesian_grid

"Analytic electric field, paired below with gridded B and force fields."
E_analytic(r) = SA[0, 5.0e-11, 0]  # [V/m]

@testset "mixed type fields" begin
    x, y, z = uniform_grid(15, 20, 25)
    B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
    F = fill(0.0, 3, length(x), length(y), length(z)) # [N]

    B[3, :, :, :] .= 1.0e-11
    F[1, :, :, :] .= 9.10938356e-42

    grid = cartesian_grid(x, y, z)

    x0 = [0.0, 0.0, 0.0] # initial position, [m]
    u0 = [0.0, 1.0, 0.0] # initial velocity, [m/s]
    stateinit = [x0..., u0...]
    tspan = (0.0, 1.0)

    param = prepare(grid, E_analytic, B, F; species = Electron)
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1, 2])
    @test isapprox(sol[1, end], 1.5324506733176166, rtol = 1.0e-2) &&
        isapprox(sol[2, end], -2.8156470047903706, rtol = 1.0e-2)
end

end # module test_mixed_fields
