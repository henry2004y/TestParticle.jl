module test_traceproblem

using Test
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using SciMLBase

@testset "TraceProblem Unified Interface" begin
    uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
    uniform_E(x) = SA[1.0e-9, 0.0, 0.0]

    x0 = [0.0, 0.0, 0.0]
    v0 = [1.0e5, 0.0, 0.0]
    stateinit_vec = [x0..., v0...]
    stateinit_sa = SA[x0..., v0...]
    tspan = (0.0, 1.0e-6)
    dt = 1.0e-9

    param = prepare(uniform_E, uniform_B; species = Proton)

    @testset "In-place solving (Vector)" begin
        prob_tp = TraceProblem(stateinit_vec, tspan, param)
        prob_ode = ODEProblem(trace!, stateinit_vec, tspan, param)

        sol_tp = solve(prob_tp, Tsit5())
        sol_ode = solve(prob_ode, Tsit5())

        @test sol_tp.t ≈ sol_ode.t
        @test sol_tp.u[end] ≈ sol_ode.u[end] rtol = 1.0e-12

        sol_vern = solve(prob_tp, Vern7())
        @test sol_vern.t[end] ≈ tspan[2]
    end

    @testset "Out-of-place solving (StaticArray)" begin
        prob_tp = TraceProblem(stateinit_sa, tspan, param)
        prob_ode = ODEProblem(trace, stateinit_sa, tspan, param)

        sol_tp = solve(prob_tp, Tsit5())
        sol_ode = solve(prob_ode, Tsit5())

        @test sol_tp.t ≈ sol_ode.t
        @test sol_tp.u[end] ≈ sol_ode.u[end] rtol = 1.0e-12
    end

    @testset "Boris solver compatibility" begin
        prob_tp = TraceProblem(stateinit_vec, tspan, param)
        sol_boris = solve(prob_tp, Boris(); dt)

        @test sol_boris.t[end] ≈ tspan[2]
        @test length(sol_boris.u) > 1
    end

    @testset "Custom ODE function" begin
        prob_rel = TraceProblem(trace_relativistic!, stateinit_vec, tspan, param)
        sol_rel = solve(prob_rel, Tsit5())
        @test sol_rel.t[end] ≈ tspan[2]
    end

    @testset "remake behavior" begin
        prob_tp = TraceProblem(stateinit_vec, tspan, param)
        @test SciMLBase.isinplace(prob_tp) == true

        # Remake with StaticArray switches to out-of-place
        prob_remade_sa = remake(prob_tp; u0 = stateinit_sa)
        @test SciMLBase.isinplace(prob_remade_sa) == false
        sol_sa = solve(prob_remade_sa, Tsit5())
        @test sol_sa.t[end] ≈ tspan[2]

        # Remake with Vector remains in-place
        prob_remade_vec = remake(prob_remade_sa; u0 = stateinit_vec)
        @test SciMLBase.isinplace(prob_remade_vec) == true
        sol_vec = solve(prob_remade_vec, Tsit5())
        @test sol_vec.t[end] ≈ tspan[2]
    end

    @testset "Ensemble solving" begin
        prob_func_test(prob, ctx) = remake(
            prob; u0 = [prob.u0[1:3]..., ctx.sim_id * 1.0e4, 0.0, 0.0]
        )
        prob_tp = TraceProblem(stateinit_vec, tspan, param; prob_func = prob_func_test)

        trajectories = 4
        sols_serial = solve(
            prob_tp, Tsit5(), EnsembleSerial(); trajectories
        )
        sols_threads = solve(
            prob_tp, Tsit5(), EnsembleThreads(); trajectories
        )

        @test length(sols_serial.u) == trajectories
        @test length(sols_threads.u) == trajectories
        for i in 1:trajectories
            @test sols_serial.u[i].u[end] ≈ sols_threads.u[i].u[end]
        end

        sols_boris = solve(
            prob_tp, Boris(), EnsembleSerial(); dt, trajectories
        )
        @test length(sols_boris.u) == trajectories
    end
end

end # module test_traceproblem
