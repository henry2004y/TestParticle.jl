if !isdefined(Main, :test_common)
    include("test_common.jl")
end

module test_traceproblem

using Test
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using SciMLBase
using ..test_common: uniform_B, uniform_Ex

@testset "TraceProblem Unified Interface" begin
    x0 = [0.0, 0.0, 0.0]
    v0 = [1.0e5, 0.0, 0.0]
    stateinit_vec = [x0..., v0...]
    stateinit_sa = SA[x0..., v0...]
    tspan = (0.0, 1.0e-6)
    dt = 1.0e-9

    param = prepare(uniform_Ex, uniform_B; species = Proton)

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
        u0_orig = copy(prob_tp.u0)

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
            @test sols_threads.u[i].prob.u0[4] ≈ i * 1.0e4
        end
        @test prob_tp.u0 == u0_orig

        sols_boris = solve(
            prob_tp, Boris(), EnsembleSerial(); dt, trajectories
        )
        @test length(sols_boris.u) == trajectories

        sols_boris_threads = solve(
            prob_tp, Boris(), EnsembleThreads(); dt, trajectories, safetycopy = false
        )
        @test length(sols_boris_threads.u) == trajectories
        for i in 1:trajectories
            @test sols_boris.u[i].u[end] ≈ sols_boris_threads.u[i].u[end]
            @test sols_boris_threads.u[i].prob.u0[4] ≈ i * 1.0e4
        end
        @test prob_tp.u0 == u0_orig
    end

    @testset "Float32 tracing" begin
        u0_f32 = Float32[x0..., v0...]
        tspan_f32 = (0.0f0, 1.0f-6)
        dt_f32 = 1.0f-9

        # Automatic promotion when passing Float32 u0 with Float64 param
        prob_f32 = TraceProblem(u0_f32, tspan, param)
        @test eltype(prob_f32.u0) === Float32
        @test eltype(prob_f32.tspan) === Float32
        @test prob_f32.p[1] isa Float32
        @test prob_f32.p[2] isa Float32

        # ODE solver (Tsit5) in Float32
        sol_ode_f32 = solve(prob_f32, Tsit5())
        @test eltype(sol_ode_f32.t) === Float32
        @test eltype(sol_ode_f32.u[end]) === Float32

        # Boris solver in Float32
        sol_boris_f32 = solve(prob_f32, Boris(); dt = dt_f32)
        @test eltype(sol_boris_f32.t) === Float32
        @test eltype(sol_boris_f32.u[end]) === Float32

        # remake switching from Float64 to Float32
        prob_f64 = TraceProblem(stateinit_vec, tspan, param)
        prob_switched = remake(prob_f64; u0 = u0_f32)
        @test eltype(prob_switched.u0) === Float32
        @test eltype(prob_switched.tspan) === Float32
        @test prob_switched.p[1] isa Float32

        # prepare with type = Float32
        param_f32 = prepare(uniform_Ex, uniform_B; species = Proton, type = Float32)
        @test param_f32[1] isa Float32
        @test param_f32[2] isa Float32
        prob_direct_f32 = TraceProblem(u0_f32, tspan_f32, param_f32)
        sol_direct = solve(prob_direct_f32, Boris(); dt = dt_f32)
        @test eltype(sol_direct.t) === Float32
        @test eltype(sol_direct.u[end]) === Float32
    end
end

end # module test_traceproblem
