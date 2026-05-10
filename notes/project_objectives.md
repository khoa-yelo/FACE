# FACE — project objectives

**Date:** 2026-05-09

The dataset-building pipeline in this repo (PhagesDB metadata + EM/plaque
images + Pharokka annotations + ESM-2 protein embeddings) is the input layer
for the FACE model itself, which is not yet in this repo. These are the
objectives the model is being built for.

## 1. Identify morphology- and plaque-determining genes (incl. novel ones)

Use joint embeddings / cross-attention between EM-image features and the
phage's per-protein ESM-2 embeddings to identify which genes drive observed
morphology and plaque phenotypes. Beyond known structural classes
(`tail`, `head and packaging`, etc.), the goal is to surface **novel
candidate determinants** — proteins that the cross-attention weights light
up on but that PHROG categorizes as `unknown function`.

Why this matters here: the recent Pharokka audit showed PHROG is
**under-labeling** real tail proteins (e.g. PHROG 6654, see
[`pharokka_tail_annotation_gap.md`](pharokka_tail_annotation_gap.md)).
ESM-2 embeddings carry the structural signal regardless of label, so a
joint model that doesn't pre-filter by category should pick those up.

**Inputs available now:**
- `data/derived_data/esm2_embeddings/<PhageName>.h5` — per-protein middle +
  last layer mean-pooled embeddings (5,980 phages, 639,493 proteins).
- `data/processed_data/em_thumb_clean/<PhageName>.png` — 200×200 RGB.
- `*_cds_categories.tsv` — per-protein PHROG category, useful as a weak
  label / prior, **not** a hard filter.

## 2. Phage classification from joint embeddings

Classify phages along multiple axes from the joint embedding:
- temperate vs lytic
- morphotype (SIPHO / MYO / PODO / TECTI / unknown)
- (other axes TBD — cluster, host, etc.)

This is partly a sanity check (does the joint embedding recover labels we
already have?) and partly a way to fill in **missing morphotype labels**
— `metadata.csv` has many phages with EM but no morphotype assignment
(e.g. Gremlin23, Noely, TicTac surfaced earlier). A trained classifier can
propose labels for those rows.

**Inputs available now:**
- `morphotype` and `temperate` columns in
  `data/processed_data/metadata.csv`.
- Cluster / host fields for downstream classification axes.

## 3. Cross-modal generation

Two directions, both predictive of how well FACE has learned the alignment:

- **Genome / proteome → EM image** — generate a plausible morphology image
  from a phage's protein set. Useful for phages with sequence but no EM, and
  as a generative diagnostic of what the model has learned about
  structure-from-sequence.
- **EM image → PHROG protein set** — given an EM, predict which PHROG
  profiles (or which protein embeddings) are likely present. Useful when
  no genome is available, and to highlight which PHROGs are
  morphology-deterministic.

## Cross-cutting notes

- Don't filter training data by `category=='tail'` etc. — the PHROG
  category column is lossy (see objective 1's note). Pass embeddings in
  raw and let the model learn the structural signal.
- Keep `metadata.csv` as the single source of truth for labels (morphotype,
  temperate, cluster). Re-derive splits from it; don't bake splits into
  preprocessing.

