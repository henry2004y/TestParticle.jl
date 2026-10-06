module test_raw_output

using Test
using TestParticle
import TestParticle as TP
using StaticArrays
using KernelAbstractions

uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
zero_E(x) = SA[0.0, 0.0, 0.0]
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
            base = TP.solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            same = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, u0 = states,
            )

            @test length(same.u) == N
            @test all(i -> same.u[i].u[end] == base.u[i].u[end], 1:N)
            @test_throws ArgumentError TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, u0 = states[1:(N - 1), :]
            )
        end

        # Fewer particles than the thread chunking threshold takes the serial path.
        let N = 4, (; prob, states, dt) = make_prob(N)
            base = TP.solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            raw = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, raw_output = true,
            )
            @test size(raw.u) == (N, 6, 2)
            @test max_state_diff(base, raw) == 0.0
        end
    end

    @testset "raw output reproduces every saved state" begin
        let N = 40, (; prob, dt) = make_prob(N)

            sols = TP.solve(prob, Boris(), CPU(); trajectories = N, dt)
            raw = TP.solve(prob, Boris(), CPU(); trajectories = N, dt, raw_output = true)
            @test size(raw.u) == (N, 6, NT + 1)
            @test raw.t == sols.u[1].t
            @test max_state_diff(sols, raw) == 0.0

            saveat = range(0.0, dt * NT, length = 11)
            sols_saveat = TP.solve(
                prob, Boris(), CPU(); trajectories = N, dt, saveat, save_everystep = false
            )
            raw_saveat = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, saveat, save_everystep = false, raw_output = true,
            )
            @test raw_saveat.t == sols_saveat.u[1].t
            @test max_state_diff(sols_saveat, raw_saveat) == 0.0

            sols_end = TP.solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            raw_end = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, raw_output = true,
            )
            @test size(raw_end.u) == (N, 6, 2)
            @test max_state_diff(sols_end, raw_end) == 0.0

            # With only the final state saved, the raw output aliases the state buffer.
            sols_final = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, save_start = false,
            )
            raw_final = TP.solve(
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

    @testset "raw output combines with u0 and particle sorting" begin
        let N = 40, (; prob, states, dt) = make_prob(N)
            reference = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false, raw_output = true,
            )

            both = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false,
                u0 = states, raw_output = true,
            )
            @test both.u == reference.u

            sorted = TP.solve(
                prob, Boris(), CPU();
                trajectories = N, dt, save_everystep = false,
                sort_particles = true, raw_output = true,
            )
            @test sorted.u == reference.u
            @test sorted.t == reference.t
        end
    end

    @testset "raw output keeps the requested element type" begin
        let N = 40, (; prob, dt) = make_prob(N; T = Float32)
            sols = TP.solve(
                prob, Boris(), CPU(); trajectories = N, dt, save_everystep = false
            )
            raw = TP.solve(
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
