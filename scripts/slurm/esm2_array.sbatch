#!/bin/bash
#SBATCH --job-name=esm2_embed
#SBATCH --partition=gpu,owners
#SBATCH --gres=gpu:1
#SBATCH --array=0-9%5
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm_logs/esm2_logs/task_%A_%a.out
#SBATCH --error=slurm_logs/esm2_logs/task_%A_%a.err

# 10 array tasks × 600 genomes/task = 6,000 capacity (covers ~5,980 phages,
# tail tasks no-op when their chunk is empty). %5 caps concurrent GPU jobs;
# with 10 tasks @ throttle 5 the array completes in roughly 2 waves.
# Each task loads the ESM-2 model once and streams the chunk through it via
# the manifest mode in esm2_embed.py — far less overhead than per-phage loads.
# Bump the throttle (number after %) up if GPU availability is good; lower it
# if you're getting evicted on the owners partition.

cd /scratch/users/khoang99/face
mkdir -p slurm_logs/esm2_logs data/derived_data/esm2_embeddings

export CHUNK_SIZE=600
export ENV_NAME=esm2
export MODEL=facebook/esm2_t33_650M_UR50D
export BATCH_SIZE=8
export MAX_LENGTH=1024

bash scripts/slurm/run_esm2_chunk.sh
