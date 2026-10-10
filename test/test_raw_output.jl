if !isdefined(Main, :test_common)
    include("test_common.jl")
end

module test_raw_output

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using KernelAbstractions
using ..test_common: uniform_B, zero_E

uniform_B32(x) = SA[0.0f0, 0.0f0, 1.0f-8]
zero_E32(x) = SA[0.0f0, 0.0f0, 0.0f0]

const DT = 1.0e-8
const NT = 200

"A distinct initial state per particle, so any reordering is visible."
function initial_states(N; T = Float64)
    states = zeros(T, N, 6)
    for i in 1:N
        states[i, 4] = i * 1.0e4
    end
    return states
end

function make_prob(N; T = Float64)
    dt = T(DT)
    states = initial_states(N; T)
    prob_func = (prob, ctx) -> remake(prob; u0 = collect(states[ctx.sim_id, :]))
    param = if T === Float32
        prepare(zero_E32, uniform_B32; species = Proton, type = Float32)
    else
        prepare(zero_E, uniform_B; species = Proton)
    end
    prob = TraceProblem(
        collect(states[1, :]), (zero(T), dt * NT), param; prob_func
    )
    return (; prob, states, dt = DT, N)
end

"The largest absolute difference between the saved states of a solution and raw output."
function max_state_diff(sols, raw)
    nout = size(raw.u, 3)
    diff = 0.0
    for i in eachindex(sols.u), j in 1:nout, k in 1:6
        diff = max(diff, abs(Float64(sols.u[i].u[j][k]) - Float64(raw.u[i, k, j])))
    end
    return diff
end

@testset "raw output" begin
    @testset "a u0 matrix replaces prob_func" begin
        let N = 40, (; prob, states, dt) = make_prob(N)
            base = solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            same = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, u0 = states,
            )

            @test length(same.u) == N
            @test all(i -> same.u[i].u[end] == base.u[i].u[end], 1:N)
            @test_throws ArgumentError solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, u0 = states[1:(N - 1), :]
            )
        end

        # Fewer particles than the thread chunking threshold takes the serial path.
        let N = 4, (; prob, states, dt) = make_prob(N)
            base = solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            raw = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, raw_output = true,
            )
            @test size(raw.u) == (N, 6, 2)
            @test max_state_diff(base, raw) == 0.0
        end
    end

    @testset "raw output reproduces every saved state" begin
        let N = 40, (; prob, dt) = make_prob(N)

            sols = solve(prob, Boris(), CPU(); trajectories = N, dt)
            raw = solve(prob, Boris(), CPU(); trajectories = N, dt, raw_output = true)
            @test size(raw.u) == (N, 6, NT + 1)
            @test raw.t == sols.u[1].t
            @test max_state_diff(sols, raw) == 0.0

            saveat = range(0.0, dt * NT, length = 11)
            sols_saveat = solve(
                prob, Boris(), CPU(); trajectories = N, dt, saveat, save_everystep = false
            )
            raw_saveat = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, saveat, save_everystep = false, raw_output = true,
            )
            @test raw_saveat.t == sols_saveat.u[1].t
            @test max_state_diff(sols_saveat, raw_saveat) == 0.0

            sols_end = solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            raw_end = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, raw_output = true,
            )
            @test size(raw_end.u) == (N, 6, 2)
            @test max_state_diff(sols_end, raw_end) == 0.0

            # With only the final state saved, the raw output aliases the state buffer.
            sols_final = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, save_start = false,
            )
            raw_final = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, save_start = false,
                raw_output = true,
            )
            @test size(raw_final.u) == (N, 6, 1)
            @test length(raw_final.t) == 1
            @test raw_final.t == sols_final.u[1].t
            @test max_state_diff(sols_final, raw_final) == 0.0
        end
    end

    @testset "raw output combines with u0" begin
        let N = 40, (; prob, states, dt) = make_prob(N)
            reference = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, raw_output = true,
            )

            both = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false,
                u0 = states, raw_output = true,
            )
            @test both.u == reference.u
            @test both.t == reference.t
        end
    end

    @testset "raw output keeps the requested element type" begin
        let N = 40, (; prob, dt) = make_prob(N; T = Float32)
            sols = solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            raw = solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, raw_output = true,
            )
            @test eltype(raw.u) === Float32
            @test eltype(raw.t) === Float32
            @test max_state_diff(sols, raw) == 0.0
        end
    end
end

end
