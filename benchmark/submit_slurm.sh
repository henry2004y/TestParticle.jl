#!/bin/bash
#SBATCH --job-name=TP_multinode
#SBATCH --constraint=cpu
#SBATCH --output=res_scaling_%j.txt
#SBATCH --error=err_scaling_%j.txt
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=00:10:00
#SBATCH --qos=debug

# One Julia worker per node, each running EnsembleThreads over the cores of that node.
# The master process stays single-threaded: SciMLBase drives EnsembleSplitThreads through
# Distributed's asynchronous `pmap`, which needs no OS threads on the master.
export JULIA_NUM_THREADS=1

# Add `--account=<project>` here or on the sbatch command line if your cluster requires it.
#
# Override --nodes, --ntasks-per-node and --cpus-per-task from the command line to sweep
# the scaling curve, e.g. `sbatch --nodes=4 submit_slurm.sh`.
echo "Nodes: $SLURM_JOB_NUM_NODES, tasks: $SLURM_NTASKS, cpus/task: $SLURM_CPUS_PER_TASK"
echo "Master threads: $JULIA_NUM_THREADS"

module load julia

# SlurmClusterManager.jl starts the workers, so the master is launched without srun.
BENCH_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
julia --project="$BENCH_DIR" "$BENCH_DIR/run_scaling_slurm.jl"
