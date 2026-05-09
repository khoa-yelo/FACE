#!/bin/bash
#SBATCH --job-name=pharokka
#SBATCH --partition=normal,owners
#SBATCH --array=0-199%50
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm_logs/pharokka_logs/task_%A_%a.out
#SBATCH --error=slurm_logs/pharokka_logs/task_%A_%a.err

# 200 array tasks × 30 genomes/task = 6,000 capacity (covers 5,980 genomes,
# tail tasks just no-op when chunk is empty). %50 caps concurrent tasks at 50.
# Adjust the throttle (number after %) up/down based on partition availability.

cd /scratch/users/khoang99/face
mkdir -p slurm_logs/pharokka_logs

export CHUNK_SIZE=30
export THREADS=4

bash scripts/slurm/run_pharokka_chunk.sh
