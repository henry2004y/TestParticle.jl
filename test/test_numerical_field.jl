module test_numerical_field

# Load the shared fixtures into Main so that `using ..test_common` resolves even
# when this file is run on its own.
if !isdefined(Main, :test_common)
    Base.include(Main, joinpath(@__DIR__, "test_common.jl"))
end

using Test
using TestParticle
import TestParticle as TP
using OrdinaryDiffEq
using StaticArrays
using Meshes: CartesianGrid, StructuredGrid
using Random
using ..test_common: uniform_grid, cartesian_grid, zero_field

"Stop the trace once the particle leaves a sphere of radius 0.8."
isoutside(u, p, t) = hypot(u[1], u[2], u[3]) > 0.8

"Rescale the initial state by a random factor, for `EnsembleProblem`."
prob_func(prob, ctx) = remake(prob; u0 = rand(ctx.rng) * prob.u0)

"Uniform field B₀ ẑ sampled on a logarithmically spaced spherical grid."
function spherical_field()
    r = logrange(0.1, 10.0, length = 4)
    θ = range(0, π, length = 8)
    ϕ = range(0, 2π, length = 8)

    B₀ = 1.0e-8 # [T]
    B = zeros(3, length(r), length(θ), length(ϕ))

    for (iθ, θ_val) in enumerate(θ)
        sinθ, cosθ = sincos(θ_val)
        B[1, :, iθ, :] .= B₀ * cosθ
        B[2, :, iθ, :] .= -B₀ * sinθ
    end
    return r, θ, ϕ, B
end

@testset "numerical field" begin
    x, y, z = uniform_grid(15, 20, 25)
    B = fill(0.0, 3, length(x), length(y), length(z)) # [T]
    E = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]

    B[3, :, :, :] .= 10.0e-9
    E[3, :, :, :] .= 5.0e-10

    grid = cartesian_grid(x, y, z)

    x0 = [0.0, 0.0, 0.0] # initial position, [m]
    u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
    stateinit = [x0..., u0...]
    tspan = (0.0, 1.0)

    param = prepare(;
        E = build_interpolator(E, (x, y, z); bc = TestParticle.FillExtrap(NaN)),
        B = build_interpolator(B, (x, y, z); bc = TestParticle.FillExtrap(NaN)),
        species = Ion(; q = 1, m = 16)
    ) # O+
    @test param[2] ≈ 16 * TP.mᵢ

    callback = TerminateOutside(isoutside)

    param = prepare(;
        E = build_interpolator(E, (x, y, z); bc = TestParticle.FillExtrap(NaN)),
        B = build_interpolator(B, (x, y, z); bc = TestParticle.FillExtrap(NaN))
    )
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(
        prob, Tsit5();
        reltol = 1.0e-8, abstol = 1.0e-8, save_idxs = [1], callback
    )
    # There are numerical differences on x86 and ARM platforms!
    @test sol[1, end] ≈ 0.79411 rtol = 1.0e-2

    param = prepare(;
        E = build_interpolator(E, (x, y, z); bc = TestParticle.FillExtrap(NaN), order = 3),
        B = build_interpolator(B, (x, y, z); bc = TestParticle.FillExtrap(NaN), order = 3)
    )
    prob = remake(prob; p = param)
    sol = solve(
        prob, Tsit5();
        reltol = 1.0e-8, abstol = 1.0e-8, save_idxs = [1], callback
    )
    @test sol[1, end] ≈ 0.79411 rtol = 1.0e-2

    # GC prepare
    stateinit_gc,
        param_gc = prepare_gc(
        stateinit;
        E = build_interpolator(E, (x, y, z); bc = TestParticle.FillExtrap(NaN), order = 1),
        B = build_interpolator(B, (x, y, z); bc = TestParticle.FillExtrap(NaN), order = 1),
        species = Proton
    )
    @test stateinit_gc[2] ≈ -1.0439684914853153
    # μ = m * v_perp^2 / (2B)
    @test param_gc[3] ≈ 8.36310961845e-20 rtol = 1.0e-6

    param = prepare(;
        E = build_interpolator(E, grid; bc = TestParticle.FillExtrap(NaN)),
        B = build_interpolator(B, grid; bc = TestParticle.FillExtrap(NaN))
    )
    prob = ODEProblem(trace!, stateinit, tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1])
    @test length(sol.u) == 8 && isapprox(sol[1, end], 0.8539409515568538, rtol = 1.0e-2)

    trajectories = 10
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sol = solve(
        ensemble_prob, Tsit5(), EnsembleThreads();
        trajectories = trajectories, save_idxs = [1], seed = 42
    )
    xs = getindex.(sol.u[10].u, 1)
    @test xs[7] ≈ 0.35505190092227473 rtol = 1.0e-6

    stateinit = SA[x0..., u0...]

    prob = ODEProblem(trace, stateinit, tspan, param)
    sol = solve(prob, Tsit5(); save_idxs = [1])
    @test length(sol.u) == 8 && isapprox(sol[1, end], 0.8539409515568538, rtol = 1.0e-2)

    # Nonuniform spherical grid
    r, θ, ϕ, B_sph = spherical_field()
    param = prepare(;
        B = build_interpolator(
            B_sph, (r, θ, ϕ);
            bc = (
                TestParticle.FillExtrap(NaN), TestParticle.FillExtrap(NaN),
                TestParticle.PeriodicBC(; endpoint = :inclusive),
            ), gridtype = StructuredGrid
        )
    )
    # Check field interpolation Bz
    @test param[4](SA[1.0, 1.0, 1.0])[3] ≈ 9.888387888463716e-9

    # Vector inputs for the spherical grid follow the same path
    param_vec = prepare(;
        B = build_interpolator(
            B_sph, (collect(r), collect(θ), collect(ϕ));
            bc = (
                TestParticle.FillExtrap(NaN), TestParticle.FillExtrap(NaN),
                TestParticle.PeriodicBC(; endpoint = :inclusive),
            ), gridtype = StructuredGrid
        ),
        species = Ion(; m = 16, q = 1)
    )
    @test param_vec[4](SA[1.0, 1.0, 1.0])[3] ≈ 9.888387888463716e-9

end

end # module test_numerical_field
