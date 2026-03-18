# SLURM EnsembleSplitThreads Scaling benchmark
#
# This script measures the performance of EnsembleSplitThreads on a SLURM cluster.
# It uses SlurmClusterManager.jl to manage workers across nodes.
#
# Usage:
#   julia --project benchmark/run_scaling_slurm.jl

using Distributed
using SlurmClusterManager

# Add workers via SLURM
# Reserve one task for the master Julia process to avoid oversubscription.
# SlurmManager() reads SLURM_NTASKS from the environment.
ntasks = parse(Int, get(ENV, "SLURM_NTASKS", "1"))
if ntasks > 1
    ENV["SLURM_NTASKS"] = string(ntasks - 1)
end

cpus_per_task = get(ENV, "SLURM_CPUS_PER_TASK", "1")
addprocs(SlurmManager(); exeflags=["--threads=$cpus_per_task"])

@everywhere begin
    using TestParticle
    using StaticArrays
    using Statistics
end

using Printf

# Benchmark settings
const N_PARTICLES = 16384
const N_SAMPLES = 5

println("Number of workers: ", nworkers())
@everywhere println("Hello from worker $(myid()) on $(gethostname()) with $(Threads.nthreads()) threads")

@everywhere begin
    # Trace Problem Setup (consistent with run_scaling_threads.jl)
    uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
    uniform_E(x) = SA[0.0, 0.0, 0.0]

    prob_func(prob, i, repeat) =
        remake(
        prob;
        u0 = [prob.u0[1], prob.u0[2], prob.u0[3], (i / 1000.0) * 1.0e5, 0.0, 0.0]
    )
end

param = prepare(uniform_E, uniform_B; species = Proton)

x0 = [0.0, 0.0, 0.0]
v0 = [1.0e5, 0.0, 0.0]
stateinit = [x0..., v0...]
tspan = (0.0, 1.0e-3)
dt = 1.0e-9

prob_multi = TraceProblem(stateinit, tspan, param; prob_func = prob_func)

println("\nStarting EnsembleSplitThreads benchmark...")
println("Particles: $N_PARTICLES, samples per point: $N_SAMPLES")

# Warmup
TestParticle.solve(
    prob_multi, EnsembleSplitThreads();
    trajectories = 10, dt, savestepinterval = 10000
)

# Timed samples
sample_times = Float64[]
for i in 1:N_SAMPLES
    t = @elapsed TestParticle.solve(
        prob_multi, EnsembleSplitThreads();
        trajectories = N_PARTICLES, dt, savestepinterval = 10000
    )
    @printf("  Sample %d: %.3f s\n", i, t)
    push!(sample_times, t)
end

t_med = median(sample_times)
@printf("\nMedian execution time: %.3f s\n", t_med)

# Save result to CSV for combined plotting later if needed
using DelimitedFiles
results_file = joinpath(@__DIR__, "slurm_scaling.csv")
writedlm(results_file, [nworkers() t_med], ',')
println("Results saved to $results_file")
