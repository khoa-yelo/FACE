#!/bin/bash
#SBATCH --job-name=phagedb_dl
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=slurm_logs/phagedb_download_%j.out
#SBATCH --error=slurm_logs/phagedb_download_%j.err

cd /scratch/users/khoang99/face
mkdir -p slurm_logs

python3 scripts/download/download_phagedb.py \
    --output-dir data/downloaded_data/phagesdb_dataset \
    --all-phages \
    --page-size 100
