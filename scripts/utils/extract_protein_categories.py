"""
Extract per-protein PHROG functional categories from pharokka results.

For each <PhageName>/ under pharokka_results/, joins phanotate.faa headers
against <PhageName>_cds_final_merged_output.tsv and writes:

    <PhageName>/<PhageName>_cds_categories.tsv
        protein_id<TAB>category

Categories are PHROG top-level classes ("tail", "head and packaging",
"unknown function", etc.) — same values pharokka emits in the `category`
column of its merged output TSV.

Usage:
    python scripts/utils/extract_protein_categories.py
    python scripts/utils/extract_protein_categories.py --phage 20ES
    python scripts/utils/extract_protein_categories.py --force
"""

import argparse
import csv
import sys
from pathlib import Path

DEFAULT_RESULTS_DIR = Path("data/derived_data/pharokka_results")


def parse_faa_ids(faa_path):
    ids = []
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def load_categories(merged_tsv):
    """Return {gene_id: category} from a pharokka *_cds_final_merged_output.tsv."""
    out = {}
    with open(merged_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            out[row["gene"]] = row["category"]
    return out


def process_phage(phage_dir, force=False):
    name = phage_dir.name
    faa = phage_dir / "phanotate.faa"
    merged = phage_dir / f"{name}_cds_final_merged_output.tsv"
    out_path = phage_dir / f"{name}_cds_categories.tsv"

    if not faa.is_file() or faa.stat().st_size == 0:
        return "skip_no_faa", 0, 0
    if not merged.is_file():
        return "skip_no_merged", 0, 0
    if out_path.is_file() and out_path.stat().st_size > 0 and not force:
        return "skip_done", 0, 0

    ids = parse_faa_ids(faa)
    cats = load_categories(merged)

    missing = 0
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["protein_id", "category"])
        for pid in ids:
            cat = cats.get(pid)
            if cat is None:
                missing += 1
                cat = "NA"
            w.writerow([pid, cat])

    return "ok", len(ids), missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR,
                    help="pharokka_results directory (default: %(default)s)")
    ap.add_argument("--phage", help="run for a single phage name only")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing _cds_categories.tsv outputs")
    args = ap.parse_args()

    if not args.results_dir.is_dir():
        print(f"ERROR: {args.results_dir} not found", file=sys.stderr)
        sys.exit(1)

    if args.phage:
        dirs = [args.results_dir / args.phage]
    else:
        dirs = sorted(p for p in args.results_dir.iterdir() if p.is_dir())

    counts = {"ok": 0, "skip_done": 0, "skip_no_faa": 0, "skip_no_merged": 0}
    total_proteins = 0
    total_missing = 0
    for d in dirs:
        status, n, missing = process_phage(d, force=args.force)
        counts[status] = counts.get(status, 0) + 1
        total_proteins += n
        total_missing += missing
        if missing:
            print(f"[warn] {d.name}: {missing}/{n} proteins had no merged-TSV row "
                  f"(written as NA)", file=sys.stderr)

    print(f"phages: ok={counts['ok']}  done={counts['skip_done']}  "
          f"no_faa={counts['skip_no_faa']}  no_merged={counts['skip_no_merged']}")
    print(f"proteins written: {total_proteins}  missing-category: {total_missing}")


if __name__ == "__main__":
    main()
