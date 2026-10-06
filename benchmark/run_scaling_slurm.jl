# Multi-node strong scaling benchmark for EnsembleSplitThreads and EnsembleDistributed.
#
# Layout: one Julia worker per node, each worker running `EnsembleThreads` over the
# cores of that node. The master process stays single-threaded on purpose: SciMLBase
# drives both ensemble algorithms through Distributed's asynchronous `pmap`, which
# needs no OS threads on the master.
#
# A single invocation is one point of the scaling curve, because the node count is
# fixed by the allocation. Sweep node counts with run_scaling_multinode.jl, or submit
# jobs by hand:
#   sbatch --nodes=4 submit_slurm.sh
#
# Usage:
#   julia --project=benchmark benchmark/run_scaling_slurm.jl           # inside a SLURM job
#   TP_LOCAL_WORKERS=2 julia --project=benchmark run_scaling_slurm.jl  # local smoke test

using Distributed
using TestParticle
using StaticArrays
using Statistics
using Printf
using DelimitedFiles

const N_PARTICLES = parse(Int, get(ENV, "TP_N_PARTICLES", "16384"))
const N_SAMPLES = parse(Int, get(ENV, "TP_N_SAMPLES", "3"))
const N_WARMUP = parse(Int, get(ENV, "TP_N_WARMUP", "16"))
const DT = 1.0e-9
const N_STEPS = parse(Int, get(ENV, "TP_N_STEPS", "1000000"))
const TSPAN = (0.0, DT * N_STEPS)
# 100 output steps per particle, matching run_scaling_threads.jl.
const SAVEAT = TSPAN[2] / 100

# Workers started without `--project` would fall back to the default environment, so
# the active project is propagated explicitly.
if haskey(ENV, "SLURM_NTASKS")
    using SlurmClusterManager
    cpus_per_task = get(ENV, "SLURM_CPUS_PER_TASK", "1")
    addprocs(
        SlurmManager();
        exeflags = ["--threads=$cpus_per_task", "--project=$(Base.active_project())"]
    )
else
    n_local = parse(Int, get(ENV, "TP_LOCAL_WORKERS", "0"))
    if n_local > 0
        addprocs(n_local; exeflags = ["--project=$(Base.active_project())"])
    end
end

# The problem and the algorithm are serialized to every worker, so everything they
# close over has to be defined there as well: `TestParticle.Field` stores the type of
# the field function, and that type is looked up by name on the receiving process.
@everywhere using TestParticle
@everywhere using StaticArrays

@everywhere uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
@everywhere uniform_E(x) = SA[0.0, 0.0, 0.0]
@everywhere prob_func(prob, ctx) = remake(
    prob;
    u0 = [prob.u0[1], prob.u0[2], prob.u0[3], (ctx.sim_id / 1000.0) * 1.0e5, 0.0, 0.0]
)

const N_NODES = parse(Int, get(ENV, "SLURM_JOB_NUM_NODES", "1"))

println("="^80)
println("Multi-node Boris ensemble scaling")
println(
    "Nodes: $N_NODES | workers: $(nworkers()) | master threads: $(Threads.nthreads())"
)
println("Particles: $N_PARTICLES, steps per particle: $N_STEPS")
for (k, w) in enumerate(workers())
    host, nthreads = fetch(@spawnat w (gethostname(), Threads.nthreads()))
    println("  worker $k (pid $w) on $host with $nthreads threads")
end
if nworkers() == 0
    @warn "No workers: both ensemble algorithms degenerate to serial execution."
end
println("="^80)

worker_threads = [fetch(@spawnat w Threads.nthreads()) for w in workers()]
threads_per_worker = isempty(worker_threads) ? 1 : maximum(worker_threads)

param = prepare(uniform_E, uniform_B; species = Proton)

stateinit = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
prob_multi = TraceProblem(stateinit, TSPAN, param; prob_func)

# Millions of particle-steps per second.
throughput(t_s) = N_PARTICLES * N_STEPS / (t_s * 1.0e6)

"""
    measure(ensemblealg, prob)

Median wall time in seconds over `N_SAMPLES` runs, after one untimed warmup that also
compiles the solver on every worker.
"""
function measure(ensemblealg, prob)
    TestParticle.solve(
        prob, Boris(), ensemblealg; trajectories = N_WARMUP, dt = DT, saveat = SAVEAT
    )
    samples = Float64[]
    for _ in 1:N_SAMPLES
        @everywhere GC.gc()
        t = @elapsed TestParticle.solve(
            prob, Boris(), ensemblealg;
            trajectories = N_PARTICLES, dt = DT, saveat = SAVEAT
        )
        push!(samples, t)
    end
    return median(samples)
end

ensemblealgs = [
    ("SplitThreads", EnsembleSplitThreads()),
    ("Distributed", EnsembleDistributed()),
]

results = NamedTuple[]
for (name, ensemblealg) in ensemblealgs
    t = measure(ensemblealg, prob_multi)
    push!(results, (; alg = name, time = t))
    @printf("  %-14s %9.3f s   %9.2f M particle-steps/s\n", name, t, throughput(t))
end

# One file per job so that concurrently running jobs never interleave their output.
job_id = get(ENV, "SLURM_JOB_ID", "local")
csv_path = joinpath(@__DIR__, "multinode_$(job_id).csv")
header = [
    "nodes" "workers" "threads_per_worker" "ensemble_alg" "particles" "median_s" "msteps_per_s"
]
open(csv_path, "w") do io
    writedlm(io, header, ',')
    for r in results
        writedlm(
            io,
            [N_NODES nworkers() threads_per_worker r.alg N_PARTICLES r.time throughput(r.time)],
            ',',
        )
    end
end
println("Results written to $csv_path")
