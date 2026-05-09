#!/bin/bash
#SBATCH --job-name=stdimg
#SBATCH --partition=normal,owners
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm_logs/standardize_images_%j.out
#SBATCH --error=slurm_logs/standardize_images_%j.err

cd /scratch/users/khoang99/face
mkdir -p slurm_logs

export PYTHONNOUSERSITE=1
micromamba run -n face python3 scripts/imaging/standardize_images.py --workers 8
