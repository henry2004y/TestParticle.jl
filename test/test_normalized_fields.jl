module test_normalized_fields

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using Test
using TestParticle
using TestParticle: ClampExtrap, WrapExtrap
import TestParticle as TP
using OrdinaryDiffEq
using StaticArrays
using ..test_common: uniform_grid, cartesian_grid

const B₀ = 10.0e-9  # [T]
const E₀ = 1.0 * B₀ # [V/m]; U₀ = 1 m/s

"A quarter gyroperiod, in normalized time."
const tspan = (0.0, 0.5π)

const x0 = [0.0, 0.0, 0.0] # initial position [l₀]
const u0 = [1.0, 0.0, 0.0] # initial velocity [v₀]
const stateinit = [x0..., u0...]

const EXPECTED = 0.9999998697180689

"Zero field in normalized units; the fields below are already divided by E₀."
zero_normalized(sz...) = zeros(3, sz...)

@testset "normalized fields" begin
    # 3D
    x, y, z = uniform_grid(15, 20, 25)
    E = zero_normalized(length(x), length(y), length(z))
    B = zero_normalized(length(x), length(y), length(z))
    B[3, :, :, :] .= 1.0

    param = prepare(;
        E = build_interpolator(E, (x, y, z); bc = TestParticle.FillExtrap(NaN)),
        B = build_interpolator(B, (x, y, z); bc = TestParticle.FillExtrap(NaN)),
        m = 1,
        q = 1
    )
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Vern9(); save_idxs = [1])
    @test length(sol.u) == 7 && isapprox(sol[1, end], 1.0, rtol = 1.0e-2)

    # 2D
    x, y, _ = uniform_grid(15, 20, 25)
    E = zero_normalized(length(x), length(y))
    B = zero_normalized(length(x), length(y))
    B[3, :, :] .= 1.0
    grid = cartesian_grid(x, y)

    # periodic BC
    param = prepare(;
        E = build_interpolator(E, grid; bc = WrapExtrap()),
        B = build_interpolator(B, grid; bc = WrapExtrap()),
        species = Proton
    )
    @test param[3] isa TP.Field
    @test startswith(repr(param[3]), "Field with")
    param = prepare(;
        E = build_interpolator(E, (x, y); bc = WrapExtrap()),
        B = build_interpolator(B, (x, y); bc = WrapExtrap()),
        m = 1,
        q = 1
    )
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1])
    @test length(sol.u) == 9 && isapprox(sol[1, end], EXPECTED, rtol = 1.0e-2)

    prob = ODEProblem(trace, SA[stateinit...], tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1])
    @test length(sol.u) == 9 && isapprox(sol[1, end], EXPECTED, rtol = 1.0e-2)

    # Because the field is uniform, the order of interpolation does not matter.
    param = prepare(;
        E = build_interpolator(E, grid; bc = TestParticle.FillExtrap(NaN), order = 3),
        B = build_interpolator(B, grid; bc = TestParticle.FillExtrap(NaN), order = 3),
        m = 1,
        q = 1
    )
    prob = remake(prob; p = param)
    sol = solve(prob, Tsit5(); save_idxs = [1])
    @test length(sol.u) == 9 && isapprox(sol[1, end], EXPECTED, rtol = 1.0e-2)

    # 1D
    x, _, _ = uniform_grid(15, 20, 25)
    E = zero_normalized(length(x))
    B = zero_normalized(length(x))
    B[3, :] .= 1.0

    param = prepare(;
        E = build_interpolator(E, (x,); bc = ClampExtrap()),
        B = build_interpolator(B, (x,); bc = ClampExtrap()),
        m = 1,
        q = 1
    )
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1])
    @test length(sol.u) == 9 && isapprox(sol[1, end], EXPECTED, rtol = 1.0e-2)
end

end # module test_normalized_fields
