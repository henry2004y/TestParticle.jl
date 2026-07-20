# # Distributed Parallelization
#
# This demo shows the multi-process ensemble parallelization strategies:
# `EnsembleDistributed()` (trajectories split across separate Julia worker
# processes) and `EnsembleSplitThreads()` (a hybrid that also uses threads
# inside each worker). Unlike the [Multithreaded parallelization](@ref) demo,
# these require starting extra worker processes and loading the needed
# packages on every worker before solving.
#
# Because the demo must run as its own Julia session with workers, the code
# below is a complete script. Run it from a terminal (e.g. `julia
# demo_parallel_distributed.jl`), not from the single-process documentation
# build.

using Distributed

# Add workers only if none exist yet, so the script is safe to re-run in the
# same session or when Julia is started with `-p`. On a cluster you would
# typically add one worker per physical core or per node.
if nworkers() == 1
    addprocs(min(4, Sys.CPU_THREADS))
end

# Every package the workers touch must be loaded on every worker.
@everywhere begin
    using TestParticle
    import TestParticle as TP
    using OrdinaryDiffEq
    using StaticArrays
end

# ## Simulation Setup
#
# The physics is identical to the multithreaded demo: electrons in a uniform
# magnetic field with Maxwellian initial velocities. All definitions that are
# evaluated on the workers live inside `@everywhere` so every process sees them.

@everywhere begin
    const Bmag = 0.01
    uniform_B(x) = SA[0.0, 0.0, Bmag]
    zero_E = ZeroField()

    param = prepare(zero_E, uniform_B; species = Electron)
    q2m = TP.get_q2m(param)

    const tperiod = 2π / (abs(q2m) * Bmag)
    const vth = 1.0e5

    tspan = (0.0, 20 * tperiod)
    dt = tperiod / 12
    stateinit = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # `prob_func` runs on the worker that owns each trajectory, so it must use
    # only worker-available functions/packages. We draw the initial velocity
    # from a 3D Maxwellian (independent Gaussian components, std = vth) using
    # the per-trajectory RNG `ctx.rng` for reproducibility.
    function prob_func(prob, ctx)
        v = SVector{3}(
            vth * randn(ctx.rng),
            vth * randn(ctx.rng),
            vth * randn(ctx.rng),
        )
        return remake(prob; u0 = vcat(SVector{3}(prob.u0[1:3]), v))
    end

    prob_ode = ODEProblem(trace, stateinit, tspan, param)
    ensemble_ode = EnsembleProblem(prob_ode; prob_func, safetycopy = false)

    prob_boris = TraceProblem(stateinit, tspan, param; prob_func)
end

trajectories = 512
seed = 1234;

# ## SciML Solver: `EnsembleDistributed`
#
# With workers available, passing `EnsembleDistributed()` distributes the
# trajectories across all worker processes via `pmap`.

sol_dist = solve(
    ensemble_ode, Vern6(), EnsembleDistributed();
    trajectories, saveat = dt, seed
);

# ## SciML Solver: `EnsembleSplitThreads` (hybrid)
#
# `EnsembleSplitThreads()` distributes trajectories across workers *and* uses
# threads inside each worker. Start Julia with `-t auto` (or set
# `JULIA_NUM_THREADS`) so each worker has multiple threads.

sol_split = solve(
    ensemble_ode, Vern6(), EnsembleSplitThreads();
    trajectories, saveat = dt, seed
);

# ## Boris Solver (same `ensemblealg` argument)
#
# The native Boris pusher accepts the identical `ensemblealg`; it only needs a
# fixed timestep `dt`.

sol_boris = TP.solve(
    prob_boris, Boris(), EnsembleDistributed();
    dt, trajectories, savestepinterval = 1, seed
);

# ## Controlling Work Granularity
#
# Distributed solvers accept `batch_size`, the number of trajectories sent to a
# worker per chunk. The default is `max(1, trajectories ÷ nworkers())`. For very
# cheap trajectories a larger `batch_size` reduces communication overhead, while
# a smaller one improves load balancing. The keyword is passed like any other
# `solve` keyword:

batch_size = 64
sol_batched = solve(
    ensemble_ode, Vern6(), EnsembleDistributed();
    trajectories, saveat = dt, seed, batch_size
);

println("Distributed trajectories:    $(length(sol_dist.u))")
println("SplitThreads trajectories:   $(length(sol_split.u))")
println("Boris (distributed) trajectories: $(length(sol_boris.u))")
println("Distributed (batch_size = $batch_size) trajectories: $(length(sol_batched.u))")

# Clean up the workers when you are done:
#
# ```julia
# rmprocs(workers())
# ```
#
# For single-session, shared-memory parallelism (no extra processes), see the
# [Multithreaded parallelization](@ref) demo.
