# FACE

Phage genome × electron-micrograph alignment.

**Status:** this repository currently contains only the **dataset-building
pipeline**. The FACE model itself (image-feature × proteome alignment) is
**in progress** and not yet part of this repo.

The pipeline builds a multimodal dataset of bacteriophages from
[PhagesDB](https://phagesdb.org): genome FASTAs, EM and plaque images, plus
uniform Pharokka annotations and ESM-2 protein embeddings. It is designed as
the input layer for the forthcoming FACE model, which will identify
structural / plaque-morphology-determining genes from morphology.

## Layout

```
.
├── README.md
├── scripts/                                # All entry points
│   ├── download/
│   │   └── download_phagedb.py             #   PhagesDB downloader (records, images, FASTAs)
│   ├── imaging/
│   │   └── standardize_images.py           #   thumbnail normalizer (200x200 RGB PNG)
│   ├── plm/
│   │   └── esm2_embed.py                   #   ESM-2 per-protein embeddings -> HDF5
│   ├── utils/
│   │   ├── compare_pharokka_vs_genbank.py  #   single-phage Pharokka vs SEA-PHAGES diff
│   │   ├── batch_compare.py                #   stratified N-phage comparison
│   │   └── convert_to_png.py               #   generic image -> PNG
│   ├── visualization/
│   │   ├── generate_morph_grids.py         #   morphotype 3x3 example grids
│   │   └── visualize_metadata.py           #   metadata distribution plots
│   └── slurm/                              #   Slurm batch entry points
│       ├── download_all_phages.sbatch
│       ├── pharokka_array.sbatch
│       ├── run_pharokka_chunk.sh
│       ├── standardize_images.sbatch
│       ├── esm2_array.sbatch
│       └── run_esm2_chunk.sh
├── data/                                   # (gitignored) all data, populated at runtime
│   ├── downloaded_data/                    #   raw data from external sources
│   │   ├── phagesdb_dataset/               #     PhagesDB downloads (images + FASTAs + metadata)
│   │   
│   ├── processed_data/                     #   cleaned/derived outputs feeding the pipeline
│   │   ├── metadata.csv                    #     canonical per-phage metadata
│   │   ├── em_thumb_clean/                 #     standardized EM thumbnails
│   │   └── plaque_thumb_clean/             #     standardized plaque thumbnails
│   └── derived_data/                       #   model outputs / intermediate artifacts
│       ├── pharokka_results/               #     per-phage Pharokka annotations
│       ├── pharokka_fasta_list.txt         #     deterministic chunking list
│       └── esm2_embeddings/                #     per-phage HDF5 protein embeddings
├── databases/                              # External DB references
│   ├── phagesdb_api.txt                    #   PhagesDB OpenAPI schema
│   └── pharokka_db/                        #   (gitignored) Pharokka DB after install
├── figures/                                # (gitignored) plots and grids
│   ├── morph_examples/                     #   morphotype example grids
│   └── phagesdb_metadata_summary.png
├── slurm_logs/                             # (gitignored) per-job Slurm logs
│   ├── pharokka_logs/
│   └── esm2_logs/
└── test/                                   # tests
```

All scripts assume the repository root as the working directory. Data paths
are relative (`data/...`, `databases/...`) so the pipeline is portable across
machines.

## Pipeline

### 1. Database setup

Download every PhagesDB record: metadata, EM micrographs, plaque images,
their thumbnails, and the genome FASTA.

```
sbatch scripts/slurm/download_all_phages.sbatch
```

Output: `data/downloaded_data/phagesdb_dataset/{em_images,plaque_images,em_thumb_images,plaque_thumb_images,fasta_files}/`
plus `metadata.csv` (one row per phage with image paths, GenBank/RefSeq
accessions, cluster, host, morphotype, etc.). The canonical metadata copy
used by downstream tools lives at `data/processed_data/metadata.csv`.

### 2. Sequence annotation (Pharokka)

Annotate every sequenced phage in parallel via a Slurm array. Default mode:
PHANOTATE + PHROG + VFDB + CARD + tRNAscan-SE + Aragorn + minced + INPHARED.

```
# one-time setup (env name `face`, database written under databases/pharokka_db/)
micromamba create -n face -c conda-forge -c bioconda pharokka python=3.11
PYTHONNOUSERSITE=1 micromamba run -n face install_databases.py -o databases/pharokka_db/

# annotate
sbatch scripts/slurm/pharokka_array.sbatch    # 200 array tasks × 30 genomes
```

`scripts/slurm/run_pharokka_chunk.sh` is the per-task worker. It is
idempotent — safe to resubmit after preemption or partial failure; genomes
whose `.gbk` already exists are skipped. Outputs land in
`data/derived_data/pharokka_results/<PhageName>/`.

### 3. Annotation evaluation (optional)

Compare Pharokka calls against the SEA-PHAGES GenBank annotation for any
phage that has a GenBank accession in `metadata.csv`.

```
# single phage
PYTHONNOUSERSITE=1 micromamba run -n face python scripts/utils/compare_pharokka_vs_genbank.py 20ES KJ410132

# stratified sample of N phages spanning N clusters
PYTHONNOUSERSITE=1 micromamba run -n face python scripts/utils/batch_compare.py -n 10
```

Reports CDS-call overlap, boundary disagreements, strand-flip events, and
function-name agreement. Caches NCBI fetches under
`data/downloaded_data/genbank_cache/`.

### 4. Image processing

Standardize EM and plaque thumbnails into a uniform tensor-ready format
(200×200 RGB PNG, alpha composited onto white, all color modes unified to
RGB, 16-bit downcast to 8-bit).

```
sbatch scripts/slurm/standardize_images.sbatch
```

Outputs land in `data/processed_data/em_thumb_clean/` and
`data/processed_data/plaque_thumb_clean/`. Unrecoverable files are logged to
`data/processed_data/standardize_failed.tsv`.

`scripts/utils/convert_to_png.py` is a generic file-by-file image-to-PNG
converter for ad-hoc use. `scripts/visualization/generate_morph_grids.py`
makes 3×3 example grids per ICTV morphotype (SIPHO/MYO/PODO/TECTI/unknown)
into `figures/morph_examples/` for sanity-checking the dataset.
`scripts/visualization/visualize_metadata.py` plots the metadata distribution
(morphotype × cluster × sequencing-status × image availability) into
`figures/phagesdb_metadata_summary.png`.

### 5. Protein embeddings (ESM-2)

For each phage with a Pharokka annotation, embed every predicted protein
(`phanotate.faa`) with ESM-2 and write per-protein mean-pooled embeddings
(middle + last transformer layers) to a per-phage HDF5 file.

```
sbatch scripts/slurm/esm2_array.sbatch    # 10 GPU array tasks × 600 genomes
```

`scripts/slurm/run_esm2_chunk.sh` is the per-task worker. Idempotent: any
phage whose `.h5` already exists is skipped. The model is loaded once per
task and reused across the chunk via the `--manifest` mode of
`scripts/plm/esm2_embed.py`. Outputs land in
`data/derived_data/esm2_embeddings/<PhageName>.h5`.

## Environment notes

All Python scripts run in the `face` micromamba environment. Pharokka brings
most dependencies (biopython, pandas, pyhmmer, mmseqs2, mash, prodigal,
phanotate, tRNAscan-SE, aragorn, minced). ESM-2 requires `transformers`,
`torch`, and `h5py` on top of that.

Always export `PYTHONNOUSERSITE=1` before invoking python in this env, to
prevent `~/.local/lib/python3.11/site-packages/` from shadowing env packages.
