"""
Compare a Pharokka annotation to the SEA-PHAGES GenBank submission for the
same phage. Report CDS-call agreement and where each annotator has features
the other doesn't.

Usage:
    python compare_pharokka_vs_genbank.py PHAGE_NAME GENBANK_ACCESSION

Requires biopython (already in the `face` env). Fetches the reference GBK
from NCBI Entrez, caches it under data/downloaded_data/genbank_cache/.
"""

import sys
import time
import argparse
from pathlib import Path
from urllib.request import Request, urlopen

from Bio import SeqIO

PHAROKKA_RESULTS = Path("data/derived_data/pharokka_results")
CACHE_DIR = Path("data/downloaded_data/genbank_cache")
EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text"
)

OVERLAP_THRESHOLD = 0.8  # IoU for calling two CDS "the same"


def fetch_genbank(accession: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out = CACHE_DIR / f"{accession}.gbk"
    if out.exists() and out.stat().st_size > 0:
        return out
    req = Request(EFETCH_URL.format(acc=accession),
                  headers={"User-Agent": "phage-compare/1.0 (academic)"})
    with urlopen(req, timeout=60) as resp:
        data = resp.read()
    out.write_bytes(data)
    time.sleep(0.34)  # NCBI rate limit
    return out


def load_cds(gbk_path: Path):
    """Return list of (start, end, strand, product, locus_tag, gene)."""
    rec = next(SeqIO.parse(gbk_path, "genbank"))
    rows = []
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        loc = feat.location
        # Use 0-based half-open coords for IoU math
        start, end = int(loc.start), int(loc.end)
        strand = "+" if loc.strand == 1 else "-"
        product = feat.qualifiers.get("product", [""])[0]
        locus = feat.qualifiers.get("locus_tag", [""])[0]
        gene = feat.qualifiers.get("gene", [""])[0]
        note = feat.qualifiers.get("note", [""])[0]
        rows.append({
            "start": start, "end": end, "strand": strand,
            "product": product, "locus_tag": locus, "gene": gene,
            "note": note,
        })
    return rec, rows


def iou(a, b):
    inter = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]))
    union = max(a["end"], b["end"]) - min(a["start"], b["start"])
    return inter / union if union else 0.0


def match_cds(pharokka, genbank):
    """For each genbank CDS, find best-IoU pharokka CDS on same strand."""
    matches = []
    used_p = set()
    for g in genbank:
        best = None
        best_iou = 0.0
        best_idx = None
        for i, p in enumerate(pharokka):
            if i in used_p:
                continue
            if p["strand"] != g["strand"]:
                continue
            score = iou(p, g)
            if score > best_iou:
                best_iou = score
                best = p
                best_idx = i
        if best_iou >= OVERLAP_THRESHOLD:
            matches.append((g, best, best_iou))
            used_p.add(best_idx)
        else:
            matches.append((g, None, best_iou))
    pharokka_only = [pharokka[i] for i in range(len(pharokka)) if i not in used_p]
    return matches, pharokka_only


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phage_name")
    parser.add_argument("accession")
    parser.add_argument("--show-mismatches", type=int, default=10,
                        help="Print up to N mismatches in each category")
    args = parser.parse_args()

    p_gbk = PHAROKKA_RESULTS / args.phage_name / f"{args.phage_name}.gbk"
    if not p_gbk.exists():
        sys.exit(f"Pharokka output not found: {p_gbk}")

    print(f"Fetching reference GenBank for {args.accession} ...")
    g_gbk = fetch_genbank(args.accession)

    p_rec, p_cds = load_cds(p_gbk)
    g_rec, g_cds = load_cds(g_gbk)

    print()
    print(f"=== {args.phage_name} ({args.accession}) ===")
    print(f"  Pharokka length: {len(p_rec.seq):>7} bp   CDS: {len(p_cds)}")
    print(f"  GenBank length:  {len(g_rec.seq):>7} bp   CDS: {len(g_cds)}")
    print()

    matches, pharokka_only = match_cds(p_cds, g_cds)
    matched = [(g, p, s) for g, p, s in matches if p is not None]
    genbank_only = [g for g, p, _ in matches if p is None]

    print(f"  Matched (IoU >= {OVERLAP_THRESHOLD}): {len(matched)}")
    print(f"  GenBank-only (Pharokka missed): {len(genbank_only)}")
    print(f"  Pharokka-only (extra calls):    {len(pharokka_only)}")
    print()

    # Function-call agreement among matched pairs
    same_function = 0
    pharokka_hypothetical = 0
    genbank_hypothetical = 0
    for g, p, _ in matched:
        gp = (g["product"] or "").lower()
        pp = (p["product"] or "").lower()
        g_hyp = ("hypothetical" in gp) or not gp
        p_hyp = ("hypothetical" in pp) or not pp
        if g_hyp:
            genbank_hypothetical += 1
        if p_hyp:
            pharokka_hypothetical += 1
        if (gp == pp) or (g_hyp and p_hyp):
            same_function += 1

    print(f"  Among {len(matched)} matched CDS:")
    print(f"    same/equivalent product:        {same_function}")
    print(f"    Pharokka hypothetical:          {pharokka_hypothetical}")
    print(f"    GenBank hypothetical:           {genbank_hypothetical}")
    print()

    n = args.show_mismatches
    if genbank_only:
        print(f"  GenBank CDS missed by Pharokka (first {min(n, len(genbank_only))}):")
        for g in genbank_only[:n]:
            print(f"    {g['locus_tag'] or '?':<14} {g['start']+1:>6}-{g['end']:<6} {g['strand']}  {g['product'] or '(no product)'}")
        print()

    if pharokka_only:
        print(f"  Pharokka CDS not in GenBank (first {min(n, len(pharokka_only))}):")
        for p in pharokka_only[:n]:
            print(f"    {p['locus_tag'] or '?':<14} {p['start']+1:>6}-{p['end']:<6} {p['strand']}  {p['product'] or '(no product)'}")
        print()


if __name__ == "__main__":
    main()
