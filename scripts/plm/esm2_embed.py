"""
ESM-2 Protein Embeddings
========================
Embed every protein in a FASTA (.faa) file with an ESM-2 model and write per-protein
mean-pooled embeddings from the middle and last transformer layers to an HDF5 file.

Usage:
    python esm2_embed.py --fasta proteins.faa --output embeddings.h5

HDF5 layout:
    /protein_id      (N,)        utf-8 strings, FASTA header IDs (text before first space)
    /middle_layer    (N, D)      float32, mean over amino acid positions of layer M
    /last_layer      (N, D)      float32, mean over amino acid positions of the final layer

Notes:
    - "Middle layer" defaults to num_hidden_layers // 2.
    - Special tokens (CLS/EOS) and padding are excluded from the mean.
    - Sequences longer than --max-length are truncated.
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer


def parse_fasta(path):
    """Yield (id, sequence) pairs from a FASTA file. ID is text before first whitespace."""
    pid, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(chunks)
                pid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if pid is not None:
            yield pid, "".join(chunks)


def batched(items, batch_size):
    batch = []
    for item in items:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


@torch.no_grad()
def embed_batch(model, tokenizer, sequences, middle_idx, last_idx, device, max_length):
    """Return (middle_means, last_means) as numpy float32 arrays of shape (B, D)."""
    enc = tokenizer(
        sequences,
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=max_length,
    )
    enc = {k: v.to(device) for k, v in enc.items()}

    out = model(**enc, output_hidden_states=True)
    hidden_states = out.hidden_states  # tuple of (B, L, D), len = num_hidden_layers + 1

    # Build a mask that keeps only real amino acid positions (drop special tokens + padding).
    attn = enc["attention_mask"]  # (B, L)
    special = torch.zeros_like(attn)
    # ESM-2 tokenizer marks CLS/EOS via special_tokens_mask in encode; recompute here.
    special_ids = set(tokenizer.all_special_ids)
    input_ids = enc["input_ids"]
    for sid in special_ids:
        special |= (input_ids == sid).long()
    aa_mask = (attn.bool() & ~special.bool()).float().unsqueeze(-1)  # (B, L, 1)
    denom = aa_mask.sum(dim=1).clamp(min=1.0)  # (B, 1)

    def pool(layer_idx):
        h = hidden_states[layer_idx]  # (B, L, D)
        return ((h * aa_mask).sum(dim=1) / denom).float().cpu().numpy()

    return pool(middle_idx), pool(last_idx)


def embed_one_file(fasta_path, out_path, model, tokenizer, middle_idx, last_idx,
                   hidden_dim, device, batch_size, max_length, model_name):
    n_total = sum(1 for _ in parse_fasta(fasta_path))
    print(f"[info] {n_total} proteins in {fasta_path}", flush=True)
    if n_total == 0:
        print(f"[warn] no proteins, skipping {fasta_path}", flush=True)
        return

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    str_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(out_path, "w") as h5:
        ids_ds = h5.create_dataset("protein_id", shape=(n_total,), dtype=str_dtype)
        mid_ds = h5.create_dataset("middle_layer", shape=(n_total, hidden_dim),
                                   dtype="float32", chunks=(min(64, n_total), hidden_dim),
                                   compression="gzip", compression_opts=4)
        last_ds = h5.create_dataset("last_layer", shape=(n_total, hidden_dim),
                                    dtype="float32", chunks=(min(64, n_total), hidden_dim),
                                    compression="gzip", compression_opts=4)
        h5.attrs["model"] = model_name
        h5.attrs["middle_layer_index"] = middle_idx
        h5.attrs["last_layer_index"] = last_idx
        h5.attrs["pooling"] = "mean over amino acid positions (special tokens excluded)"
        h5.attrs["max_length"] = max_length

        cursor = 0
        for batch in batched(parse_fasta(fasta_path), batch_size):
            ids = [pid for pid, _ in batch]
            seqs = [seq for _, seq in batch]
            mids, lasts = embed_batch(model, tokenizer, seqs, middle_idx, last_idx,
                                      device, max_length)
            n = len(ids)
            ids_ds[cursor:cursor + n] = np.array(ids, dtype=object)
            mid_ds[cursor:cursor + n] = mids
            last_ds[cursor:cursor + n] = lasts
            cursor += n

    print(f"[done] wrote {out_path}", flush=True)


def main():
    p = argparse.ArgumentParser(description="Compute ESM-2 mean-pooled protein embeddings into HDF5.")
    p.add_argument("--fasta", help="Input .faa / .fasta file (single-file mode).")
    p.add_argument("--output", help="Output .h5 file (single-file mode).")
    p.add_argument("--manifest", help="TSV file with columns: fasta_path<TAB>output_h5_path. "
                                      "Model is loaded once and reused across all entries.")
    p.add_argument("--skip-existing", action="store_true",
                   help="In manifest mode, skip entries whose output already exists with nonzero size.")
    p.add_argument("--model", default="facebook/esm2_t33_650M_UR50D",
                   help="HuggingFace model name (default: facebook/esm2_t33_650M_UR50D).")
    p.add_argument("--middle-layer", type=int, default=None,
                   help="Index of the 'middle' layer in hidden_states (0 = embeddings, "
                        "1..N = transformer blocks). Default: num_hidden_layers // 2.")
    p.add_argument("--batch-size", type=int, default=4)
    p.add_argument("--max-length", type=int, default=1024,
                   help="Max tokens per sequence (longer sequences are truncated).")
    p.add_argument("--device", default=None, help="cuda / cpu (default: auto).")
    args = p.parse_args()

    if not args.manifest and not (args.fasta and args.output):
        sys.exit("Provide either --manifest, or both --fasta and --output.")

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[info] device={device}  model={args.model}", flush=True)

    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = AutoModel.from_pretrained(args.model)
    model.eval().to(device)

    num_layers = model.config.num_hidden_layers  # transformer blocks; hidden_states has +1
    last_idx = num_layers
    middle_idx = args.middle_layer if args.middle_layer is not None else num_layers // 2
    if not (0 <= middle_idx <= num_layers):
        sys.exit(f"--middle-layer {middle_idx} out of range [0, {num_layers}]")
    hidden_dim = model.config.hidden_size
    print(f"[info] layers={num_layers}  middle_idx={middle_idx}  last_idx={last_idx}  "
          f"hidden_dim={hidden_dim}", flush=True)

    if args.manifest:
        with open(args.manifest) as f:
            jobs = [line.rstrip("\n").split("\t") for line in f if line.strip()]
        print(f"[info] manifest has {len(jobs)} entries", flush=True)
        for i, row in enumerate(jobs, 1):
            if len(row) != 2:
                print(f"[warn] manifest line {i} malformed, skipping: {row}", flush=True)
                continue
            fasta_path, out_path = row
            if args.skip_existing and Path(out_path).exists() and Path(out_path).stat().st_size > 0:
                print(f"[skip] {out_path} already exists ({i}/{len(jobs)})", flush=True)
                continue
            print(f"[run]  {i}/{len(jobs)}  {fasta_path} -> {out_path}", flush=True)
            try:
                embed_one_file(fasta_path, out_path, model, tokenizer, middle_idx, last_idx,
                               hidden_dim, device, args.batch_size, args.max_length, args.model)
            except Exception as e:
                print(f"[fail] {fasta_path}: {e}", flush=True)
    else:
        embed_one_file(args.fasta, args.output, model, tokenizer, middle_idx, last_idx,
                       hidden_dim, device, args.batch_size, args.max_length, args.model)


if __name__ == "__main__":
    main()
