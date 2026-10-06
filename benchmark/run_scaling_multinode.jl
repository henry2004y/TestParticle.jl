# Multi-node scaling sweep driver
#
# Submits one SLURM job per node count, waits for them to finish, and merges the
# per-job CSVs written by run_scaling_slurm.jl into multinode_scaling.csv. Run it from
# a login node; the measurements themselves never run there.
#
# Usage:
#   julia --project=benchmark benchmark/run_scaling_multinode.jl [nodes...]
#   julia --project=benchmark benchmark/run_scaling_multinode.jl 1 2 4 8
#
# Node counts default to DEFAULT_NODES. Extra sbatch arguments can be injected with
# TP_SBATCH_ARGS, e.g.
#   TP_SBATCH_ARGS="--qos=regular --time=00:30:00 --account=m1234" julia ...
#
# Plot the collected curve with plot_scaling_multinode.jl.

using DelimitedFiles
using Printf

const DEFAULT_NODES = [1, 2, 4, 8]
const BENCH_DIR = @__DIR__
const POLL_INTERVAL = 15
const ACTIVE_STATES = Set(
    [
        "PENDING", "RUNNING", "CONFIGURING", "SUSPENDED", "COMPLETING", "REQUEUE", "RESIZING",
    ]
)

const HEADER = [
"nodes" "workers" "threads_per_worker" "ensemble_alg" "particles" "median_s" "msteps_per_s"
]

function clear_previous_results()
    for f in readdir(BENCH_DIR)
        if startswith(f, "multinode_") && endswith(f, ".csv") && f != "multinode_scaling.csv"
            rm(joinpath(BENCH_DIR, f); force = true)
        end
    end
    return nothing
end

function submit_job(nodes, extra_args)
    script = joinpath(BENCH_DIR, "submit_slurm.sh")
    cmd = `sbatch --parsable --nodes=$nodes $extra_args $script`
    out = try
        read(cmd, String)
    catch err
        @warn "sbatch failed for $nodes node(s)" exception = err
        return nothing
    end
    m = match(r"\d+", out)
    if m === nothing
        @warn "Could not parse a job id from sbatch output" output = out
        return nothing
    end
    return m.match
end

function active_states(ids)
    out = try
        read(`squeue -j $(join(ids, ",")) -h -o %T`, String)
    catch
        return String[]
    end
    return strip.(split(strip(out), '\n'; keepempty = false))
end

function wait_for_jobs(ids)
    while !isempty(ids)
        states = active_states(ids)
        any(in(ACTIVE_STATES), states) || return nothing
        sleep(POLL_INTERVAL)
    end
    return nothing
end

"""
    collect_results() -> Vector{Vector{Any}}

Read every per-job `multinode_*.csv` and return the rows sorted by node count and
ensemble algorithm.
"""
function collect_results()
    rows = Vector{Vector{Any}}()
    for f in sort(readdir(BENCH_DIR))
        startswith(f, "multinode_") && endswith(f, ".csv") || continue
        f == "multinode_scaling.csv" && continue
        data = readdlm(joinpath(BENCH_DIR, f), ',')
        for i in 2:size(data, 1)
            push!(rows, Any[data[i, j] for j in 1:size(data, 2)])
        end
    end
    return sort!(rows; by = r -> (Float64(r[1]), String(r[4])))
end

function report(rows)
    baseline = Dict{String, Float64}()
    for r in rows
        alg = String(r[4])
        haskey(baseline, alg) || (baseline[alg] = Float64(r[6]))
    end

    println()
    @printf(
        "%-6s | %-8s | %-13s | %-10s | %-12s | %-8s\n",
        "Nodes", "Workers", "Ensemble", "Time (s)", "Msteps/s", "Speedup"
    )
    println("-"^72)
    for r in rows
        alg = String(r[4])
        @printf(
            "%-6d | %-8d | %-13s | %10.3f | %12.2f | %7.2fx\n",
            Int(Float64(r[1])), Int(Float64(r[2])), alg,
            Float64(r[6]), Float64(r[7]), baseline[alg] / Float64(r[6])
        )
    end

    out_path = joinpath(BENCH_DIR, "multinode_scaling.csv")
    open(out_path, "w") do io
        writedlm(io, [HEADER "speedup"], ',')
        for r in rows
            alg = String(r[4])
            writedlm(
                io,
                [r[1] r[2] r[3] alg r[5] r[6] r[7] baseline[alg] / Float64(r[6])],
                ',',
            )
        end
    end
    println("\nSaved merged results to $out_path")
    return nothing
end

function main()
    node_counts = isempty(ARGS) ? DEFAULT_NODES : parse.(Int, ARGS)
    extra_args = split(get(ENV, "TP_SBATCH_ARGS", ""))

    clear_previous_results()

    ids = String[]
    for nodes in node_counts
        id = submit_job(nodes, extra_args)
        id === nothing && continue
        push!(ids, id)
        println("Submitted $nodes node(s) as job $id")
    end
    isempty(ids) && error("No jobs were submitted; see the warnings above.")

    wait_for_jobs(ids)

    rows = collect_results()
    isempty(rows) &&
        error("No results collected; check res_scaling_*.txt and err_scaling_*.txt.")
    report(rows)
    println("Run plot_scaling_multinode.jl to plot the curve.")
    return nothing
end

if isempty(PROGRAM_FILE) || abspath(PROGRAM_FILE) == @__FILE__
    main()
end
