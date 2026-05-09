"""Batch-compare Pharokka annotations vs SEA-PHAGES GenBank for N phages.

Picks N phages stratified across clusters, fetches GenBanks (cached),
and tabulates per-phage and pooled stats:
    matched / gb_only_function / gb_only_hyp / pharokka_only_function /
    pharokka_only_hyp / strand_flipped (gb call overlapped by opposite-strand
    pharokka call) / boundary_disagreement (function-only gb call has
    same-strand pharokka overlap but IoU < 0.8).
"""
import csv
import sys
import argparse
import random
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_pharokka_vs_genbank import (
    load_cds, fetch_genbank, PHAROKKA_RESULTS, OVERLAP_THRESHOLD
)

METADATA_CSV = Path('data/processed_data/metadata.csv')


def is_hyp(prod):
    p = (prod or "").lower()
    return ("hypothetical" in p) or not p


def overlap_bp(a, b):
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]))


def classify(p_cds, g_cds, threshold=OVERLAP_THRESHOLD):
    """Return per-category counts for one phage."""
    used_p = set()
    matched = 0
    same_function_matched = 0

    # Pass 1: strict same-strand IoU>=threshold matching
    for g in g_cds:
        best_idx = None
        best_iou = 0.0
        for i, p in enumerate(p_cds):
            if i in used_p or p["strand"] != g["strand"]:
                continue
            inter = overlap_bp(p, g)
            union = max(p["end"], g["end"]) - min(p["start"], g["start"])
            iou = inter / union if union else 0.0
            if iou > best_iou:
                best_iou = iou
                best_idx = i
        if best_iou >= threshold and best_idx is not None:
            used_p.add(best_idx)
            matched += 1
            gp = (g["product"] or "").lower()
            pp = (p_cds[best_idx]["product"] or "").lower()
            if (gp == pp) or (is_hyp(g["product"]) and is_hyp(p_cds[best_idx]["product"])):
                same_function_matched += 1
            g["_matched"] = True
        else:
            g["_matched"] = False

    # Classify remaining GB-only calls
    gb_only_func_truly_missed = 0
    gb_only_func_boundary = 0
    gb_only_func_strandflip = 0
    gb_only_hyp = 0
    for g in g_cds:
        if g["_matched"]:
            continue
        if is_hyp(g["product"]):
            gb_only_hyp += 1
            continue
        # Has function — check overlap with any pharokka CDS
        best_same = 0
        best_opp = 0
        for p in p_cds:
            ov = overlap_bp(p, g)
            if p["strand"] == g["strand"]:
                best_same = max(best_same, ov)
            else:
                best_opp = max(best_opp, ov)
        glen = g["end"] - g["start"]
        if best_same / glen > 0.05:
            gb_only_func_boundary += 1
        elif best_opp / glen > 0.5:
            gb_only_func_strandflip += 1
        else:
            gb_only_func_truly_missed += 1

    # Pharokka-only
    pharokka_only_func = 0
    pharokka_only_hyp = 0
    for i, p in enumerate(p_cds):
        if i in used_p:
            continue
        if is_hyp(p["product"]):
            pharokka_only_hyp += 1
        else:
            pharokka_only_func += 1

    return {
        "p_cds": len(p_cds),
        "g_cds": len(g_cds),
        "matched": matched,
        "same_function": same_function_matched,
        "gb_only_func_missed": gb_only_func_truly_missed,
        "gb_only_func_boundary": gb_only_func_boundary,
        "gb_only_func_strandflip": gb_only_func_strandflip,
        "gb_only_hyp": gb_only_hyp,
        "pharokka_only_func": pharokka_only_func,
        "pharokka_only_hyp": pharokka_only_hyp,
    }


def pick_phages(n, seed=42):
    """Pick n phages stratified across clusters: pharokka done + has GB acc."""
    csv.field_size_limit(10**8)
    by_cluster = defaultdict(list)
    with open(METADATA_CSV, newline='') as f:
        for row in csv.DictReader(f):
            if (row.get('seq_finished') != 'True'
                or not (row.get('genbank_accession') or '').strip()):
                continue
            name = row['phage_name']
            gbk = PHAROKKA_RESULTS / name / f"{name}.gbk"
            if not gbk.exists() or gbk.stat().st_size == 0:
                continue
            cluster = (row.get('cluster') or 'Unknown').strip() or 'Unknown'
            by_cluster[cluster].append((name, row['genbank_accession'].strip()))
    rng = random.Random(seed)
    clusters = sorted(by_cluster, key=lambda c: -len(by_cluster[c]))
    picked = []
    while len(picked) < n and clusters:
        for c in list(clusters):
            if not by_cluster[c]:
                clusters.remove(c)
                continue
            choice = rng.choice(by_cluster[c])
            by_cluster[c].remove(choice)
            picked.append((c, *choice))
            if len(picked) >= n:
                break
    return picked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    picks = pick_phages(args.n, seed=args.seed)
    print(f"Comparing {len(picks)} phages (one per cluster, seed={args.seed}):\n")

    rows = []
    pooled = defaultdict(int)
    for cluster, name, acc in picks:
        try:
            p_rec, p_cds = load_cds(PHAROKKA_RESULTS / name / f"{name}.gbk")
            g_rec, g_cds = load_cds(fetch_genbank(acc))
            stats = classify(p_cds, g_cds)
            stats.update(name=name, acc=acc, cluster=cluster,
                         length=len(p_rec.seq))
            rows.append(stats)
            for k, v in stats.items():
                if isinstance(v, int) and k not in ("length",):
                    pooled[k] += v
        except Exception as e:
            print(f"  [{name}/{acc}] error: {e}")

    # Per-phage table
    hdr = ["cluster", "phage", "acc", "len", "P", "GB", "match", "=fn",
           "GBmiss", "GBbound", "GBstrand", "GBhyp", "Pextra", "Phyp"]
    print(" ".join(f"{h:<8}" for h in hdr))
    for r in rows:
        print(f"{r['cluster']:<8} {r['name']:<8} {r['acc']:<8} "
              f"{r['length']:<8} {r['p_cds']:<8} {r['g_cds']:<8} "
              f"{r['matched']:<8} {r['same_function']:<8} "
              f"{r['gb_only_func_missed']:<8} {r['gb_only_func_boundary']:<8} "
              f"{r['gb_only_func_strandflip']:<8} {r['gb_only_hyp']:<8} "
              f"{r['pharokka_only_func']:<8} {r['pharokka_only_hyp']:<8}")

    # Pooled summary
    print()
    print("Pooled across all phages:")
    total_p = pooled['p_cds']
    total_g = pooled['g_cds']
    print(f"  Pharokka CDS total:                 {total_p}")
    print(f"  GenBank  CDS total:                 {total_g}")
    print(f"  Matched (IoU≥{OVERLAP_THRESHOLD}, same strand):    {pooled['matched']} "
          f"({100*pooled['matched']/total_g:.1f}% of GB)")
    print(f"  Of matched, same product:           {pooled['same_function']} "
          f"({100*pooled['same_function']/pooled['matched']:.1f}% of matched)")
    print()
    print(f"  GenBank-only (Pharokka misses):")
    print(f"    truly missed (no overlap):        {pooled['gb_only_func_missed']}")
    print(f"    boundary disagreement:            {pooled['gb_only_func_boundary']}")
    print(f"    strand-flipped overlap:           {pooled['gb_only_func_strandflip']}")
    print(f"    hypothetical (uninformative):     {pooled['gb_only_hyp']}")
    print()
    print(f"  Pharokka-only (extra calls):")
    print(f"    with function (PHROG hit):        {pooled['pharokka_only_func']}")
    print(f"    hypothetical:                     {pooled['pharokka_only_hyp']}")


if __name__ == "__main__":
    main()
