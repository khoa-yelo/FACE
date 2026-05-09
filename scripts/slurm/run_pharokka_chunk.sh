#!/bin/bash
# Annotate one chunk of phage FASTAs with Pharokka.
# Called per-task by pharokka_array.sbatch — uses SLURM_ARRAY_TASK_ID
# to pick which slice of the FASTA list this task processes.
#
# Layout:
#   FASTA_DIR/PhageName.fasta  -> RESULTS_DIR/PhageName/PhageName.gbk (etc.)
# A task is idempotent: if RESULTS_DIR/PhageName/PhageName.gbk exists
# and has nonzero size, that genome is skipped.

set -euo pipefail

FASTA_DIR="/scratch/users/khoang99/face/data/downloaded_data/phagesdb_dataset/fasta_files"
RESULTS_DIR="/scratch/users/khoang99/face/data/derived_data/pharokka_results"
DB_DIR="/scratch/users/khoang99/face/databases/pharokka_db"
LIST_FILE="/scratch/users/khoang99/face/data/derived_data/pharokka_fasta_list.txt"
CHUNK_SIZE="${CHUNK_SIZE:-30}"
THREADS="${THREADS:-4}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set" >&2
    exit 1
fi

mkdir -p "$RESULTS_DIR"

# Build the FASTA list once if missing (sorted for deterministic chunking)
if [[ ! -s "$LIST_FILE" ]]; then
    ls "$FASTA_DIR"/*.fasta | sort > "$LIST_FILE"
fi

start=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE + 1 ))
end=$(( start + CHUNK_SIZE - 1 ))
chunk=$(sed -n "${start},${end}p" "$LIST_FILE")

if [[ -z "$chunk" ]]; then
    echo "[task $SLURM_ARRAY_TASK_ID] empty chunk (lines $start-$end), exiting"
    exit 0
fi

echo "[task $SLURM_ARRAY_TASK_ID] processing lines $start-$end ($(echo "$chunk" | wc -l) genomes)"

# Avoid leaking ~/.local site-packages into pharokka's python
export PYTHONNOUSERSITE=1

success=0
skipped=0
failed=0
while IFS= read -r fasta; do
    [[ -z "$fasta" ]] && continue
    name=$(basename "$fasta" .fasta)
    out="$RESULTS_DIR/$name"
    gbk="$out/${name}.gbk"

    if [[ -s "$gbk" ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    rm -rf "$out"
    if micromamba run -n face pharokka.py \
        -i "$fasta" \
        -o "$out" \
        -d "$DB_DIR" \
        -t "$THREADS" \
        -p "$name" \
        -f \
        > "$out.log" 2>&1; then
        success=$((success + 1))
        rm -f "$out.log"
    else
        failed=$((failed + 1))
        echo "[task $SLURM_ARRAY_TASK_ID] FAIL: $name (see $out.log)"
    fi
done <<< "$chunk"

echo "[task $SLURM_ARRAY_TASK_ID] done: $success succeeded, $skipped skipped, $failed failed"
exit $(( failed > 0 ? 1 : 0 ))
