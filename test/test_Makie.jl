# test Makie recipes

module TestModule
using TestParticle, OrdinaryDiffEq, StaticArrays
using CairoMakie
using Test

# a basic case for TestParticle
B(xu, t) = SA[0, 0, 1e-8*sin(20π*t)+1e-8]
E(xu) = SA[0.0, 1e-7, 0.0]
x0 = [10.0, 10.0, 0.0] # initial position, [m]
u0 = [1000.0, 0.0, 10.0] # initial velocity, [m/s]
tspan = (0.0, 0.1)
stateinit = [x0..., u0...]
param = prepare(E, B; species=Electron)

sol = let prob = ODEProblem(trace!, stateinit, tspan, param)
    solve(prob, Vern6())
end

sol_boris = let dt = 0.01
    prob = TraceProblem(stateinit, tspan, param)
    TestParticle.solve(prob; dt)[1]
end

@testset "orbit" begin
    fig = orbit(sol)
    @test fig isa Figure

    fig = orbit(sol, idxs=(1, 2))
    @test fig isa Figure

    f(xu) = xu[1]
    fig = orbit(sol, idxs=f, tspan=(0.0, 0.05))
    @test fig isa Figure

    fig = orbit(sol, idxs=(1, [2, 3]))
    @test fig isa Figure

    fig = orbit(sol, idxs=([1, 2], 3))
    @test fig isa Figure

    fig = orbit(sol, idxs=([1, 2], [3, 4]))
    @test fig isa Figure

    fig = orbit(sol, idxs=[(1, 2, 3), (1, 2)], to_3d=true)
    @test fig isa Figure

    @test_throws ArgumentError orbit(sol, idxs=([1, 2], "x"))

    @test_throws ArgumentError orbit(sol, idxs=(1, "x"))
end

@testset "monitor" begin
    fig = monitor(sol)
    @test fig isa Figure
end

end