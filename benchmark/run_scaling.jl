# Parallel scaling benchmark workflow for TestParticle.
#
# One entry point covers the three parallelism axes of the package:
#
#   threads      one node, `EnsembleThreads`, sweeps Julia thread counts
#   distributed  one node, `EnsembleDistributed`, sweeps worker process counts
#   multinode    one worker per node, `EnsembleSplitThreads` against
#                `EnsembleDistributed`, swept through the SLURM queue
#
# Every point is measured by the same `point` subcommand in a fresh process, so the
# three modes differ only in how resources are provisioned, and all of them share one
# CSV schema (see ScalingCommon.jl). Measurement never plots; render figures afterwards
# with plot_scaling.jl.
#
# Usage:
#   julia --project=benchmark benchmark/run_scaling.jl threads [counts...]
#   julia --project=benchmark benchmark/run_scaling.jl distributed [counts...]
#   julia --project=benchmark benchmark/run_scaling.jl multinode [counts...]
#   julia --project=benchmark benchmark/run_scaling.jl collect
#   julia -t 8 --project=benchmark benchmark/run_scaling.jl point threads 8 out.csv
#
# Every mode is also submittable from a login node through submit_slurm.sh, which passes
# its arguments straight through:
#   sbatch submit_slurm.sh threads
#   sbatch submit_slurm.sh distributed 1 2 4 8 16 32 64 128
#   sbatch --nodes=4 submit_slurm.sh multinode
#
# `threads` and `distributed` run on a single node and write threads_scaling.csv and
# distributed_scaling.csv. `multinode` measures one node count and writes
# multinode_<jobid>.csv, one file per job. Inside an allocation it detects `SLURM_JOB_ID`
# and measures; on a login node it submits one submit_slurm.sh job per node count, waits
# for the queue to drain, and merges. To drive the submission through SLURM yourself,
# sbatch the jobs with `--nodes=N` and then run `collect` to merge the per-job CSVs.
#
# Counts default to 1, 2, 4, 8 for `multinode` and to the powers of two up to the node's
# logical core count for the other modes. Problem size and repetition come from the
# environment, see ScalingCommon.jl. Extra sbatch flags go in TP_SBATCH_ARGS, e.g.
#   TP_SBATCH_ARGS="--account=m1234 --qos=regular" julia ... run_scaling.jl multinode

using Distributed
using Printf
using TestParticle

include(joinpath(@__DIR__, "ScalingCommon.jl"))
using .ScalingCommon

const BENCH_DIR = @__DIR__
const SCRIPT_PATH = @__FILE__
const DEFAULT_NODES = [1, 2, 4, 8]
const POLL_INTERVAL = 15
const ACTIVE_STATES = Set(
    [
        "PENDING", "RUNNING", "CONFIGURING", "SUSPENDED", "COMPLETING", "REQUEUE",
        "RESIZING",
    ]
)

"""
    default_local_counts() -> Vector{Int}

Powers of two up to the node's logical core count, with that count appended when it is
not itself a power of two.
"""
function default_local_counts()
    max_count = Sys.CPU_THREADS
    counts = [2^i for i in 0:floor(Int, log2(max_count))]
    max_count in counts || push!(counts, max_count)
    return counts
end

function print_points(points)
    isempty(points) && return nothing
    speedup = speedups(points)
    println()
    @printf(
        "%-12s | %-13s | %-9s | %-10s | %-12s | %-8s\n",
        "Mode", "Ensemble", "Count", "Time (s)", "Msteps/s", "Speedup"
    )
    println("-"^76)
    for (p, s) in zip(points, speedup)
        @printf(
            "%-12s | %-13s | %-9d | %10.3f | %12.2f | %7.2fx\n",
            p.mode, p.ensemble_alg, p.count, p.median_s, throughput(p.median_s), s
        )
    end
    return nothing
end

"""
    print_environment()

Log the Julia version and the TestParticle build being measured, so a job log always
records which package the numbers belong to.
"""
function print_environment()
    progress("Julia $(VERSION), project $(something(Base.active_project(), "none"))")
    progress("TestParticle $(pkgversion(TestParticle)) at $(pathof(TestParticle))")
    return nothing
end

"""
    setup_workers(pids)

Load ScalingCommon.jl on every worker. `TestParticle` serializes `prob_func` and the
field functions by name, so the module they are defined in has to exist on each process
that receives the problem.
"""
function setup_workers(pids)
    path = joinpath(BENCH_DIR, "ScalingCommon.jl")
    @everywhere pids include($path)
    return nothing
end

function point_threads(count)
    actual = Threads.nthreads()
    actual == count || @warn "Julia started with $actual threads" requested = count
    t = measure(EnsembleThreads(), make_problem())
    return [
        ScalingPoint(
            "threads", "threads", actual, 1, actual, "Threads", N_PARTICLES, N_STEPS, t,
        ),
    ]
end

function point_distributed(count)
    pids = addprocs(count; exeflags = ["--project=$(Base.active_project())"])
    try
        setup_workers(pids)
        t = measure(EnsembleDistributed(), make_problem())
        return [
            ScalingPoint(
                "distributed", "workers", nworkers(), nworkers(), 1, "Distributed",
                N_PARTICLES, N_STEPS, t,
            ),
        ]
    finally
        rmprocs(pids)
    end
end

function point_multinode(count)
    haskey(ENV, "SLURM_JOB_NUM_NODES") ||
        error("multinode measurement must run inside a SLURM allocation")
    n_nodes = parse(Int, ENV["SLURM_JOB_NUM_NODES"])
    n_nodes == count || @warn "Allocation has $n_nodes node(s)" requested = count

    isdefined(Main, :SlurmClusterManager) || @eval using SlurmClusterManager
    cpus_per_task = get(ENV, "SLURM_CPUS_PER_TASK", "1")
    pids = addprocs(
        SlurmManager();
        exeflags = ["--threads=$cpus_per_task", "--project=$(Base.active_project())"]
    )
    try
        setup_workers(pids)
        threads_per_worker = maximum(fetch(@spawnat w Threads.nthreads()) for w in pids)
        prob = make_problem()
        results = [
            ("SplitThreads", measure(EnsembleSplitThreads(), prob)),
            ("Distributed", measure(EnsembleDistributed(), prob)),
        ]
        return [
            ScalingPoint(
                "multinode", "nodes", n_nodes, length(pids), threads_per_worker,
                name, N_PARTICLES, N_STEPS, t,
            ) for (name, t) in results
        ]
    finally
        rmprocs(pids)
    end
end

function point_mode(mode, count, out_path)
    progress(
        "$mode point: count=$count, particles=$N_PARTICLES, steps=$N_STEPS, " *
            "samples=$N_SAMPLES"
    )
    points = if mode == "threads"
        point_threads(count)
    elseif mode == "distributed"
        point_distributed(count)
    elseif mode == "multinode"
        point_multinode(count)
    else
        error("Unknown measurement mode: $mode")
    end
    write_points(out_path, points)
    print_points(points)
    println("Wrote $(length(points)) row(s) to $out_path")
    return nothing
end

const LOCAL_OUTPUTS = Dict(
    "threads" => "threads_scaling.csv", "distributed" => "distributed_scaling.csv"
)

function run_local_sweep(mode, counts)
    out_path = joinpath(BENCH_DIR, LOCAL_OUTPUTS[mode])
    progress("="^72)
    progress("$(uppercasefirst(mode)) sweep over counts $counts")
    progress("$N_PARTICLES particles, $N_STEPS steps, $N_SAMPLES sample(s) per point")
    progress("Each point runs in its own process; the smallest count is the slowest.")
    progress("Output: $out_path")
    progress("="^72)
    tmp = mktempdir(BENCH_DIR)
    points = ScalingPoint[]
    try
        for count in counts
            child_csv = joinpath(tmp, "$(mode)_$(count).csv")
            project = Base.active_project()
            julia = mode == "threads" ? `$(Base.julia_cmd()) -t $count` : Base.julia_cmd()
            cmd = `$julia --project=$project $SCRIPT_PATH point $mode $count $child_csv`
            progress("Measuring $(mode) point at count = $count")
            run(cmd)
            append!(points, read_points(child_csv))
        end
    finally
        rm(tmp; recursive = true, force = true)
    end
    sort!(points; by = p -> (p.count, p.ensemble_alg))
    write_points(out_path, points)
    print_points(points)
    println("\nSaved $out_path")
    return out_path
end

function clear_previous_results()
    for f in readdir(BENCH_DIR)
        startswith(f, "multinode_") && endswith(f, ".csv") || continue
        f == "multinode_scaling.csv" && continue
        rm(joinpath(BENCH_DIR, f); force = true)
    end
    return nothing
end

function submit_job(nodes, extra_args)
    script = joinpath(BENCH_DIR, "submit_slurm.sh")
    # Slurm runs the batch script from its own spool directory, so pass the script location
    # explicitly instead of relying on the directory the job was submitted from.
    location = "--export=ALL,TP_BENCH_DIR=$BENCH_DIR"
    out = try
        read(
            `sbatch --parsable --nodes=$nodes $extra_args $location $script multinode`,
            String,
        )
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
    while any(in(ACTIVE_STATES), active_states(ids))
        sleep(POLL_INTERVAL)
    end
    return nothing
end

function collect_results()
    points = ScalingPoint[]
    for f in sort(readdir(BENCH_DIR))
        startswith(f, "multinode_") && endswith(f, ".csv") || continue
        f == "multinode_scaling.csv" && continue
        append!(points, read_points(joinpath(BENCH_DIR, f)))
    end
    return points
end

"""
    collect_multinode() -> String

Merge every per-job `multinode_*.csv` into multinode_scaling.csv and print the speedup
table. This is the second half of the `multinode` driver and the whole of the `collect`
mode, which is what you run after sbatching the jobs through SLURM yourself. Duplicate
points are reported, because they mean a stale per-job CSV is still lying around.
"""
function collect_multinode()
    points = collect_results()
    isempty(points) && error(
        "No multinode_*.csv found; check res_scaling_*.txt and err_scaling_*.txt."
    )
    sort!(points; by = p -> (p.count, p.ensemble_alg))
    seen = Set{Tuple{Int, String}}()
    for p in points
        key = (p.count, p.ensemble_alg)
        key in seen &&
            @warn "Duplicate point; remove the stale multinode_*.csv" point = key
        push!(seen, key)
    end
    out_path = joinpath(BENCH_DIR, "multinode_scaling.csv")
    write_points(out_path, points)
    print_points(points)
    println("\nSaved $out_path")
    println("Render the curve with plot_scaling.jl multinode.")
    return out_path
end

function run_multinode_sweep(node_counts)
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
    collect_multinode()
    return nothing
end

function main(args = ARGS)
    print_environment()
    isempty(args) && error(
        "Usage: run_scaling.jl <threads|distributed|multinode|collect> [counts...]"
    )
    mode, rest = args[1], args[2:end]

    if mode == "point"
        length(rest) == 3 || error("Usage: run_scaling.jl point <mode> <count> <outfile>")
        point_mode(rest[1], parse(Int, rest[2]), rest[3])
    elseif haskey(LOCAL_OUTPUTS, mode)
        counts = isempty(rest) ? default_local_counts() : parse.(Int, rest)
        run_local_sweep(mode, counts)
    elseif mode == "collect"
        collect_multinode()
    elseif mode == "multinode"
        if haskey(ENV, "SLURM_JOB_ID")
            default_nodes = isempty(rest) ? "1" : rest[1]
            n_nodes = parse(Int, get(ENV, "SLURM_JOB_NUM_NODES", default_nodes))
            out = joinpath(BENCH_DIR, "multinode_$(ENV["SLURM_JOB_ID"]).csv")
            point_mode("multinode", n_nodes, out)
        else
            counts = isempty(rest) ? DEFAULT_NODES : parse.(Int, rest)
            run_multinode_sweep(counts)
        end
    else
        error("Unknown mode `$mode`; choose threads, distributed, multinode or collect.")
    end
    return nothing
end

if isempty(PROGRAM_FILE) || abspath(PROGRAM_FILE) == @__FILE__
    main()
end
