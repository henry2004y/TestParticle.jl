#!/bin/bash
#SBATCH --job-name=TP_scaling
#SBATCH --constraint=cpu
#SBATCH --output=res_scaling_%j.txt
#SBATCH --error=err_scaling_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=01:00:00
#SBATCH --qos=regular

# Batch wrapper for the three scaling modes. The first argument is the mode, the rest
# are the counts to sweep:
#
#   sbatch submit_slurm.sh threads
#   sbatch submit_slurm.sh distributed 1 2 4 8 16 32 64 128
#   sbatch --nodes=4 submit_slurm.sh multinode
#
# `threads` and `distributed` run their whole sweep inside this single job, so they want
# one node. `multinode` measures one node count, so pass `--nodes=N` and keep
# --ntasks-per-node at 1, which is what yields one worker per node. Submit the node counts
# you want and then merge them with `julia run_scaling.jl collect`, or let the driver do
# both by running `run_scaling.jl multinode 1 2 4 8` on a login node.
#
# The defaults size a single node job: the one hour covers a full sweep at the default
# problem size (about 35 minutes). Add `--account=<project>` if your cluster requires
# it, and override --nodes, --cpus-per-task, --time or --qos on the sbatch command line
# as needed.
#
# Set TP_PROJECT=<path> in the environment to choose the Julia environment that provides
# TestParticle, e.g. `TP_PROJECT=~/envs/bench sbatch submit_slurm.sh threads`.
#
# The scaling scripts are found through TP_BENCH_DIR first, then SLURM_SUBMIT_DIR, and
# finally this file's own directory. Because Slurm runs a batch script from its spool
# directory, `sbatch /somewhere/else/submit_slurm.sh multinode` needs TP_BENCH_DIR to point
# at the directory holding run_scaling.jl and ScalingCommon.jl.
#
# Pass explicit counts when you want a specific set. Otherwise the sweep goes up to the
# physical core count of the node it lands on, detected from /proc/cpuinfo, so it stops
# before hyperthreading. That detection follows the node rather than the allocation, and
# TP_MAX_THREADS overrides it.

# One Julia worker per node, each running EnsembleThreads over the cores of that node.
# The master process stays single-threaded: SciMLBase partitions the trajectories over
# workers() and only uses threads inside each worker, so the master needs no OS threads.
export JULIA_NUM_THREADS=1

if [ $# -eq 0 ]; then
    echo "Usage: sbatch submit_slurm.sh <threads|distributed|multinode> [counts...]" >&2
    exit 2
fi

echo "Nodes: $SLURM_JOB_NUM_NODES, tasks: $SLURM_NTASKS, cpus/task: $SLURM_CPUS_PER_TASK"
echo "Command: run_scaling.jl $*"

module load julia
export PATH=/global/homes/h/hyzhou/.juliaup/bin/:$PATH

# SlurmClusterManager.jl starts the workers in multinode mode, so the master is launched
# without srun.
#
# Slurm copies a batch script into its own spool directory before running it, so $0 and
# BASH_SOURCE point there rather than at the directory you submitted from. Take the script
# directory from TP_BENCH_DIR, then from the submission directory, and only fall back to
# this file's own location when running outside Slurm. run_scaling.jl's driver sets
# TP_BENCH_DIR, so it no longer depends on where you launch it from.
if [ -n "$TP_BENCH_DIR" ]; then
    BENCH_DIR="$TP_BENCH_DIR"
elif [ -n "$SLURM_SUBMIT_DIR" ]; then
    BENCH_DIR="$SLURM_SUBMIT_DIR"
else
    BENCH_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi

if [ ! -f "$BENCH_DIR/run_scaling.jl" ]; then
    echo "BENCH_DIR=$BENCH_DIR does not contain run_scaling.jl" >&2
    echo "Submit from the directory holding the scaling scripts, or set" >&2
    echo "TP_BENCH_DIR=<that directory> in the environment." >&2
    exit 2
fi

# Julia environment providing TestParticle and SlurmClusterManager. A standalone copy of
# these scripts works without a project, because Julia falls back to the default
# environment through its load path. Set TP_PROJECT=/path/to/env to pin one explicitly,
# which is what you want when the default environment holds a different TestParticle
# than the one you mean to measure.
PROJECT="${TP_PROJECT:-}"
if [ -z "$PROJECT" ] && [ -f "$BENCH_DIR/Project.toml" ]; then
    PROJECT="$BENCH_DIR"
fi

if [ -n "$PROJECT" ]; then
    echo "Project: $PROJECT"
    julia --project="$PROJECT" "$BENCH_DIR/run_scaling.jl" "$@"
else
    echo "Project: <julia default environment>"
    julia "$BENCH_DIR/run_scaling.jl" "$@"
fi
