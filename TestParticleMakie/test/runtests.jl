using TestParticleMakie, TestParticle, OrdinaryDiffEq, StaticArrays
using GLMakie
using Test

# a basic case for TestParticle
B(xu, t) = SA[0, 0, 1e-8*sin(20π*t)+1e-8]
E(xu) = SA[0.0, 1e-7, 0.0]
x0 = [10.0, 10.0, 0.0] # initial position, [m]
u0 = [1000.0, 0.0, 10.0] # initial velocity, [m/s]
tspan = (0.0, 0.1)
stateinit = [x0..., u0...]
param = prepare(E, B; species=Electron)

prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern6(); save_idxs=[1,2,3,4,5,6])

@testset "basic recipe" begin
    fig, ax, plt = plot(sol)
    @test fig isa Figure 

    fig, ax, plt = lines(sol, vars=(3, 4), to_3d=true)
    @test plt isa Lines

    fig, ax, plt = lines(sol, vars=[1, 2], to_3d=true)
    @test plt isa Lines

    fig, ax, plt = lines(sol, vars=[(1, 2), 3])
    @test plt isa Lines

    @test_throws ArgumentError lines(sol, vars="x")

    @test_throws ArgumentError lines(sol, vars=["x", "y"])

    @test_throws ArgumentError lines(sol, vars=8)
end

@testset "orbit" begin
    fig = orbit(sol)
    @test fig isa Figure

    fig = orbit(sol, vars=(1, 2))
    @test fig isa Figure

    f(xu) = xu[1]
    fig = orbit(sol, vars=f, tspan=(0.0, 0.05))
    @test fig isa Figure

    fig = orbit(sol, vars=(1, [2, 3]))
    @test fig isa Figure

    fig = orbit(sol, vars=([1, 2], 3))
    @test fig isa Figure

    fig = orbit(sol, vars=([1, 2], [3, 4]))
    @test fig isa Figure

    fig = orbit(sol, vars=[(1, 2, 3), (1, 2)], to_3d=true)
    @test fig isa Figure

    @test_throws ArgumentError orbit(sol, vars=([1, 2], "x"))

    @test_throws ArgumentError orbit(sol, vars=(1, "x"))
end

@testset "monitor" begin
    fig = monitor(sol)
    @test fig isa Figure
end