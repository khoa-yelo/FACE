#!/bin/bash
# Embed phanotate.faa proteins from each pharokka result with ESM-2.
# Called per-task by esm2_array.sbatch — uses SLURM_ARRAY_TASK_ID to pick
# which slice of the phage list this task processes.
#
# Layout:
#   RESULTS_DIR/PhageName/phanotate.faa  ->  EMBED_DIR/PhageName.h5
# A task is idempotent: outputs that already exist with nonzero size are skipped.
# The model is loaded once per task and reused across the whole chunk.

set -euo pipefail

RESULTS_DIR="/scratch/users/khoang99/face/data/derived_data/pharokka_results"
EMBED_DIR="/scratch/users/khoang99/face/data/derived_data/esm2_embeddings"
LIST_FILE="/scratch/users/khoang99/face/data/derived_data/esm2_phage_list.txt"
SCRIPT="/scratch/users/khoang99/face/scripts/plm/esm2_embed.py"
CHUNK_SIZE="${CHUNK_SIZE:-200}"
ENV_NAME="${ENV_NAME:-face}"
MODEL="${MODEL:-facebook/esm2_t33_650M_UR50D}"
BATCH_SIZE="${BATCH_SIZE:-8}"
MAX_LENGTH="${MAX_LENGTH:-1024}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set" >&2
    exit 1
fi

mkdir -p "$EMBED_DIR"

# Build the phage list once if missing (sorted for deterministic chunking).
# Only include directories whose phanotate.faa exists and is nonempty.
if [[ ! -s "$LIST_FILE" ]]; then
    find "$RESULTS_DIR" -mindepth 2 -maxdepth 2 -name phanotate.faa -size +0 \
        | sort > "$LIST_FILE"
fi

start=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE + 1 ))
end=$(( start + CHUNK_SIZE - 1 ))
chunk=$(sed -n "${start},${end}p" "$LIST_FILE")

if [[ -z "$chunk" ]]; then
    echo "[task $SLURM_ARRAY_TASK_ID] empty chunk (lines $start-$end), exiting"
    exit 0
fi

n_in_chunk=$(echo "$chunk" | wc -l)
echo "[task $SLURM_ARRAY_TASK_ID] processing lines $start-$end ($n_in_chunk genomes)"

# Build a TSV manifest (fasta\toutput.h5) for this chunk, dropping ones already done.
manifest="$EMBED_DIR/.manifest_task_${SLURM_ARRAY_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID}.tsv"
: > "$manifest"
queued=0
skipped=0
while IFS= read -r faa; do
    [[ -z "$faa" ]] && continue
    name=$(basename "$(dirname "$faa")")
    out="$EMBED_DIR/$name.h5"
    if [[ -s "$out" ]]; then
        skipped=$((skipped + 1))
        continue
    fi
    printf '%s\t%s\n' "$faa" "$out" >> "$manifest"
    queued=$((queued + 1))
done <<< "$chunk"

echo "[task $SLURM_ARRAY_TASK_ID] queued=$queued  already-done=$skipped"

if (( queued == 0 )); then
    rm -f "$manifest"
    exit 0
fi

export PYTHONNOUSERSITE=1

micromamba run -n "$ENV_NAME" python "$SCRIPT" \
    --manifest "$manifest" \
    --skip-existing \
    --model "$MODEL" \
    --batch-size "$BATCH_SIZE" \
    --max-length "$MAX_LENGTH"

status=$?
rm -f "$manifest"
echo "[task $SLURM_ARRAY_TASK_ID] python exit=$status"
exit $status
