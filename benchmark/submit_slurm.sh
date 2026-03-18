#!/bin/bash
#SBATCH --job-name=TP_scaling
#SBATCH --constraint=cpu
#SBATCH --output=res_scaling_%j.txt
#SBATCH --error=err_scaling_%j.txt
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=128
#SBATCH --time=00:10:00
#SBATCH --qos=debug

# Keep the master Julia process lightweight to allow 2 full workers on Node 1
export JULIA_NUM_THREADS=1

echo "Starting SLURM job with $SLURM_JOB_NUM_NODES nodes and $SLURM_NTASKS tasks total."
echo "Threads per task: $JULIA_NUM_THREADS"

module load julia

# Run the benchmark
# Note: We don't use srun here because SlurmClusterManager.jl handles worker creation.
# We just launch one Julia process as the master.
julia run_scaling_slurm.jl
