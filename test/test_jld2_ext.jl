module test_jld2_ext

using Test
using TestParticle
using BorisPushers
using StaticArrays
using JLD2

@testset "JLD2 Extension" begin
    tmp_dir = mktempdir()

    x0 = SA[0.0, 0.0, 0.0]
    v0 = SA[0.0, 1.0e5, 0.0]
    stateinit = [x0..., v0...]
    tspan = (0.0, 3.0e-8)
    dt = 3.0e-11

    param = prepare(; B = (x, t) -> SA[0.0, 0.0, 0.01], species = Electron)
    prob = TraceProblem(stateinit, tspan, param)
    sol = solve(prob, Boris(); dt = dt, saveat = 10 * dt)

    @testset "Single solution" begin
        file1 = joinpath(tmp_dir, "test1.h5")
        file2 = joinpath(tmp_dir, "test2.jld2")

        save_trajectory(file1, sol)
        @test isfile(file1)
        jldopen(file1) do f
            @test f["t"] ≈ sol.t
            @test f["x"] ≈ [u[1] for u in sol.u]
            @test f["y"] ≈ [u[2] for u in sol.u]
            @test f["z"] ≈ [u[3] for u in sol.u]
            @test f["vx"] ≈ [u[4] for u in sol.u]
            @test f["vy"] ≈ [u[5] for u in sol.u]
            @test f["vz"] ≈ [u[6] for u in sol.u]
            @test size(f["states"]) == (6, length(sol.t))
            @test f["states"][1, :] ≈ f["x"]
        end

        # Test reversed argument ordering
        save_trajectory(sol, file2; states_matrix = false)
        @test isfile(file2)
        jldopen(file2) do f
            @test f["t"] ≈ sol.t
            @test !haskey(f, "states")
        end
    end

    @testset "Raw vectors" begin
        file3 = joinpath(tmp_dir, "test_raw_vec.h5")
        t_raw = [0.0, 1.0, 2.0]
        u_raw = [SA[1.0, 2.0, 3.0, 4.0, 5.0, 6.0] for _ in 1:3]
        save_trajectory(file3, t_raw, u_raw)
        @test isfile(file3)
        jldopen(file3) do f
            @test f["t"] == t_raw
            @test f["x"] == [1.0, 1.0, 1.0]
            @test size(f["states"]) == (6, 3)
        end
    end

    @testset "Ensemble solution" begin
        ensemble_prob = TraceProblem(
            stateinit, tspan, param;
            prob_func = (p, ctx) -> remake(p; u0 = [p.u0[1:4]..., ctx.sim_id * 1.0e5, p.u0[6]])
        )
        sols = solve(
            ensemble_prob, Boris(), EnsembleSerial();
            dt = dt, saveat = 10 * dt, trajectories = 2
        )
        file4 = joinpath(tmp_dir, "test_ensemble.h5")
        save_trajectory(file4, sols)
        @test isfile(file4)
        jldopen(file4) do f
            @test f["num_trajectories"] == 2
            @test haskey(f, "trajectory_1")
            @test haskey(f, "trajectory_2")
            traj1 = f["trajectory_1"]
            @test traj1["t"] ≈ sols.u[1].t
            @test size(traj1["states"]) == (6, length(sols.u[1].t))
        end
    end

    @testset "Raw batch NamedTuple" begin
        file5 = joinpath(tmp_dir, "test_raw_batch.h5")
        raw_u = zeros(Float64, 4, 6, 10)
        raw_t = collect(range(0.0, 1.0; length = 10))
        save_trajectory(file5, (; u = raw_u, t = raw_t))
        @test isfile(file5)
        jldopen(file5) do f
            @test f["num_trajectories"] == 4
            @test size(f["states"]) == (4, 6, 10)
            @test f["t"] == raw_t
        end
    end

    rm(tmp_dir; recursive = true, force = true)
end

end # module test_jld2_ext
