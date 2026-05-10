# Pharokka tail-annotation gap (PHROG curation, not pharokka)

**Date:** 2026-05-09
**Tested on:** 5,980 pharokka annotations under `data/derived_data/pharokka_results/`

## Question

How many phages have no structural-gene calls from pharokka, and is that
biology (genuinely tail-less morphologies) or an annotation gap?

## Per-class structural coverage

Across all 5,980 pharokka annotations:

| Missing class | Phages | % |
|---|---|---|
| `head and packaging` | 7 | 0.1% |
| `tail` | 19 | 0.3% |
| `connector` | 694 | 11.6% |
| **all three** | **0** | **0%** |

Connector is noisy because portal/head-tail-connector proteins are commonly
folded under `head and packaging` instead — not a real gap.

## The 26 problem phages

**Missing `head and packaging` (7):**
Araxxi, DoTi, Moyashi, P1.1, P101A, P104A, PA6

**Missing `tail` (19):**
Biozilla, CrunchyBoi, Dennebes, ERHickory, Fabian, FlowerPower, GalaxyEater,
Geostin, Gremlin23, HitchHiker, Noely, Oatly, Pabst, PigHeaded,
PineapplePluto, RetrieverFever, Rideau, TicTac, Vorvolakos.

EM availability: 2/7 head-missing have EM, 15/19 tail-missing have EM
(GalaxyEater and TicTac only as raw images, no clean thumbnail).
Grids written to `figures/missing_head_em.png` and
`figures/missing_tail_em.png`.

## Morphotype distribution

| Set | Morphotypes |
|---|---|
| Missing head (with EM) | Araxxi: PODO,  Moyashi: SIPHO |
| Missing tail (with EM) | 12× PODO,  3× unlabeled (Gremlin23, Noely, TicTac) |

Initial read: PODO bias suggests genuinely short/absent tails — biology, not
a bug. **This was wrong.** See next section.

## GenBank check refutes the biology hypothesis

Of the 17 EM-having problem phages, 14 have GenBank accessions. Every one of
them has tail-related CDS in the SEA-PHAGES GenBank submission:

| phage | GB CDS-with-tail-keyword | notable hits |
|---|---|---|
| Moyashi | 14 | portal, head-to-tail adaptor & connector, **major tail protein**, minor tail protein |
| Noely | 7 | portal, head-to-tail adaptor, **tail terminator**, **major tail protein** |
| Araxxi | 3 | portal + minor tail protein |
| 11 podo phages (CrunchyBoi, Dennebes, Fabian, FlowerPower, Geostin, Gremlin23, Pabst, PineapplePluto, RetrieverFever, TicTac, Vorvolakos) | 3 each | portal + 2× minor tail protein |

So even the "obvious podos" do have minor tail proteins in GenBank — they're
just being missed by pharokka.

## Why `--genbank` mode doesn't fix it

Hypothesis: maybe PHANOTATE missed the gene boundaries. Test: re-run pharokka
on FlowerPower with `--genbank -g genbank` (use GenBank's CDS coordinates
directly).

```
pharokka.py -i data/downloaded_data/genbank_cache/MH155868.gbk \
    -o data/derived_data/pharokka_genbank_test/FlowerPower \
    -d databases/pharokka_db -t 4 -p FlowerPower \
    --genbank -g genbank -f
```

**Result: still zero `tail` calls.** Category counts barely changed (only
~10 fewer "unknown function" CDS, from PHANOTATE's extra short ORFs that
GenBank doesn't have).

**The two minor-tail-protein CDS still get matched to a PHROG profile** —
both hit `phrog=6654` — but pharokka tags them `category=unknown function`
because **PHROG 6654's curator-assigned category in the PHROG metadata is
`unknown function`, not `tail`.**

```
=== Pharokka calls at GenBank's tail coords (12172-13474, 13474-14767) ===
NNCOTDNY_CDS_0011  12173-13474  cat=unknown function  annot=hypothetical protein  phrog=6654
NNCOTDNY_CDS_0012  13475-14767  cat=unknown function  annot=hypothetical protein  phrog=6654
```

The detection signal is fine. The category label is wrong upstream.

## Implications for FACE

1. **Don't bulk re-annotate with `--genbank`.** It won't recover tail calls;
   the gap is in PHROG curation, not in pharokka's prediction.
2. **Filtering proteins by `category=='tail'` will silently lose real tail
   genes for ≥14 phages here, and likely many more across the corpus.**
   Worth running `scripts/utils/batch_compare.py` at higher N to estimate
   the corpus-wide rate.
3. **A PHROG-ID-level filter is more reliable than the category column.** If
   we curate a list of "tail-like" PHROG IDs (e.g. add 6654 manually, or
   pull from a more aggressive tail-protein resource like VOGs or the PHROG
   annotation table directly), we recover proteins the category column
   misses.
4. **ESM-2 embeddings are unaffected** — they're computed from the protein
   sequence regardless of PHROG label, so any FACE model that learns from
   embeddings + image isn't hurt by the category gap. Only models that
   pre-filter by category would be.

## Spot-check phages worth a manual look

- **Moyashi** — SIPHO morphology, no head/packaging called. Likely real
  annotation gap.
- **Noely** — unlabeled morphotype but EM clearly shows a long tail.
  Almost certainly a SIPHO; GenBank confirms (`major tail protein`,
  `tail terminator`).
- **Gremlin23, TicTac** — unlabeled morphotype, useful to assign from EM
  even independent of the tail-call question.
