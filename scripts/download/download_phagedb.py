"""
PhagesDB Image & Metadata Downloader
=====================================
Downloads EM micrographs, plaque images, and paired metadata from PhagesDB API.
Builds a structured dataset suitable for ML training.

Usage:
    python phagesdb_image_downloader.py [--output-dir OUTPUT_DIR] [--cluster CLUSTER] [--max-pages MAX]

Output structure:
    output_dir/
    ├── metadata.csv           # All phage metadata with image paths
    ├── em_images/             # Electron microscopy micrographs
    │   ├── PhageName.jpg
    │   └── ...
    ├── plaque_images/         # Plaque photos
    │   ├── PhageName.jpg
    │   └── ...
    ├── em_thumb_images/       # EM thumbnails (cropped, consistent JPGs)
    │   ├── PhageName.jpg
    │   └── ...
    ├── plaque_thumb_images/   # Plaque thumbnails (cropped, consistent JPGs)
    │   ├── PhageName.jpg
    │   └── ...
    └── fasta_files/           # Genome FASTA sequences (sequenced phages)
        ├── PhageName.fasta
        └── ...
"""

import os
import csv
import json
import time
import argparse
import urllib.request
import urllib.error
from pathlib import Path


BASE_URL = "https://phagesdb.org"
API_BASE = f"{BASE_URL}/api"

# Fields to extract for ML metadata
METADATA_FIELDS = [
    "phage_name", "morphotype", "genome_length", "gcpercent",
    "end_type", "num_ORFs", "num_tRNAs", "num_tmRNAs",
    "isolation_temperature", "found_year",
    "genbank_accession", "refseq_accession",
    "seq_finished",
    "em_file", "plaque_file", "em_thumb_file", "plaque_thumb_file",
    "fasta_file",
    # These will be extracted from nested objects
    "cluster", "subcluster", "temperate",
    "host_genus", "host_species",
    # Local file paths (filled during download)
    "em_local_path", "plaque_local_path",
    "em_thumb_local_path", "plaque_thumb_local_path",
    "fasta_local_path",
]


def fetch_json(url, retries=3, delay=1.0):
    """Fetch JSON from URL with retry logic and rate limiting."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "Accept": "application/json",
                "User-Agent": "PhagesDB-Dataset-Downloader/1.0 (academic research)",
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            print(f"  [WARN] Attempt {attempt+1}/{retries} failed for {url}: {e}")
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
    return None


def download_image(url, dest_path, retries=3):
    """Download an image file. Returns True if successful."""
    if not url:
        return False
    # Ensure full URL
    if url.startswith("/"):
        url = BASE_URL + url
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "PhagesDB-Dataset-Downloader/1.0 (academic research)",
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                with open(dest_path, "wb") as f:
                    f.write(resp.read())
            time.sleep(0.1)  # Small delay between image downloads
            return True
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(1)
    return False


def extract_metadata(phage):
    """Extract flat metadata dict from nested API response."""
    row = {}

    # Direct fields
    for field in ["phage_name", "morphotype", "genome_length", "gcpercent",
                  "end_type", "num_ORFs", "num_tRNAs", "num_tmRNAs",
                  "isolation_temperature", "found_year",
                  "genbank_accession", "refseq_accession",
                  "seq_finished",
                  "em_file", "plaque_file", "em_thumb_file", "plaque_thumb_file",
                  "fasta_file"]:
        row[field] = phage.get(field, "")

    # Nested: cluster / subcluster
    pcluster = phage.get("pcluster") or {}
    psubcluster = phage.get("psubcluster") or {}
    row["cluster"] = pcluster.get("cluster", "")
    row["temperate"] = pcluster.get("temperate", "")
    row["subcluster"] = psubcluster.get("subcluster", "")

    # Nested: host
    host_genus = phage.get("host_genus") or {}
    host_species = phage.get("host_species") or {}
    row["host_genus"] = host_genus.get("genus_name", "")
    row["host_species"] = host_species.get("species_name", "")

    # Placeholders for local paths
    row["em_local_path"] = ""
    row["plaque_local_path"] = ""
    row["em_thumb_local_path"] = ""
    row["plaque_thumb_local_path"] = ""
    row["fasta_local_path"] = ""

    return row


def fetch_all_phages(endpoint="sequenced_phages", page_size=100, max_pages=None,
                     cluster=None, filter_cluster=None):
    """
    Paginate through the PhagesDB API and yield phage records.

    If cluster is specified and filter_cluster is None, use the cluster API endpoint
    (sequenced only). If filter_cluster is set, use the main endpoint and filter
    client-side (includes unsequenced).
    """
    if cluster and not filter_cluster:
        url = f"{API_BASE}/clusters/{cluster}/phagelist/?page=1&page_size={page_size}"
    else:
        url = f"{API_BASE}/{endpoint}/?page=1&page_size={page_size}"

    page = 0
    consecutive_failures = 0
    while url:
        page += 1
        if max_pages and page > max_pages:
            print(f"  Reached max pages ({max_pages}), stopping.")
            break

        print(f"  Fetching page {page}: {url}")
        data = fetch_json(url)
        if not data:
            consecutive_failures += 1
            print(f"  [ERROR] Failed to fetch page {page} ({consecutive_failures} consecutive failures)")
            if consecutive_failures >= 5:
                print("  [ERROR] 5 consecutive page failures, stopping.")
                break
            # Try to construct the next page URL manually to skip this page
            try:
                from urllib.parse import urlparse, parse_qs, urlencode, urlunparse
                parsed = urlparse(url)
                params = parse_qs(parsed.query)
                current_page = int(params.get("page", [page])[0])
                params["page"] = [str(current_page + 1)]
                url = urlunparse(parsed._replace(query=urlencode(params, doseq=True)))
            except Exception:
                print("  [ERROR] Could not construct next page URL, stopping.")
                break
            time.sleep(2)
            continue

        consecutive_failures = 0
        results = data.get("results", [])
        print(f"  Got {len(results)} phages (total: {data.get('count', '?')})")

        for phage in results:
            if filter_cluster:
                pc = phage.get("pcluster") or {}
                if pc.get("cluster") != filter_cluster:
                    continue
            yield phage

        url = data.get("next")
        time.sleep(0.5)  # Be polite to the server


def _write_summary(summary_path, stats, cluster_filter):
    summary = {
        "total_phages": stats["total"],
        "phages_with_em_image": stats["has_em"],
        "phages_with_plaque_image": stats["has_plaque"],
        "phages_with_em_thumb": stats["has_em_thumb"],
        "phages_with_plaque_thumb": stats["has_plaque_thumb"],
        "phages_with_fasta": stats["has_fasta"],
        "em_images_downloaded": stats["em_downloaded"],
        "plaque_images_downloaded": stats["plaque_downloaded"],
        "em_thumb_downloaded": stats["em_thumb_downloaded"],
        "plaque_thumb_downloaded": stats["plaque_thumb_downloaded"],
        "fasta_downloaded": stats["fasta_downloaded"],
        "morphotype_distribution": stats["morphotypes"],
        "cluster_filter": cluster_filter,
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Download PhagesDB images + metadata")
    parser.add_argument("--output-dir", default="phagesdb_dataset",
                        help="Output directory (default: phagesdb_dataset)")
    parser.add_argument("--cluster", default=None,
                        help="Download only a specific cluster (e.g. 'A', 'K')")
    parser.add_argument("--all-phages", action="store_true",
                        help="Use /api/phages/ to include unsequenced phages (default: sequenced only)")
    parser.add_argument("--max-pages", type=int, default=None,
                        help="Max API pages to fetch (for testing)")
    parser.add_argument("--skip-images", action="store_true",
                        help="Only download metadata, skip image files")
    parser.add_argument("--page-size", type=int, default=100,
                        help="Results per API page (default: 100)")
    args = parser.parse_args()

    # Create output directories
    out_dir = Path(args.output_dir)
    em_dir = out_dir / "em_images"
    plaque_dir = out_dir / "plaque_images"
    em_thumb_dir = out_dir / "em_thumb_images"
    plaque_thumb_dir = out_dir / "plaque_thumb_images"
    fasta_dir = out_dir / "fasta_files"
    out_dir.mkdir(parents=True, exist_ok=True)
    em_dir.mkdir(exist_ok=True)
    plaque_dir.mkdir(exist_ok=True)
    em_thumb_dir.mkdir(exist_ok=True)
    plaque_thumb_dir.mkdir(exist_ok=True)
    fasta_dir.mkdir(exist_ok=True)

    stats = {"total": 0, "has_em": 0, "has_plaque": 0,
             "has_em_thumb": 0, "has_plaque_thumb": 0, "has_fasta": 0,
             "em_downloaded": 0, "plaque_downloaded": 0,
             "em_thumb_downloaded": 0, "plaque_thumb_downloaded": 0,
             "fasta_downloaded": 0,
             "morphotypes": {}, "seen_names": set()}

    print("=" * 60)
    print("PhagesDB Image & Metadata Downloader")
    print("=" * 60)
    if args.cluster:
        print(f"Cluster filter: {args.cluster}")
    print(f"Output directory: {out_dir}")
    print()

    # Incremental CSV: write header, then append each row
    csv_path = out_dir / "metadata.csv"
    csv_file = open(csv_path, "w", newline="")
    csv_writer = csv.DictWriter(csv_file, fieldnames=METADATA_FIELDS)
    csv_writer.writeheader()

    endpoint = "phages" if args.all_phages else "sequenced_phages"
    # When --all-phages + --cluster, filter client-side to include unsequenced
    filter_cluster = args.cluster if (args.all_phages and args.cluster) else None
    cluster_arg = args.cluster if not filter_cluster else None
    for phage in fetch_all_phages(endpoint=endpoint,
                                   page_size=args.page_size,
                                   max_pages=args.max_pages,
                                   cluster=cluster_arg,
                                   filter_cluster=filter_cluster):
        row = extract_metadata(phage)
        name = row["phage_name"]

        # Skip duplicate phage names
        if name in stats["seen_names"]:
            print(f"    [SKIP] Duplicate: {name}")
            continue
        stats["seen_names"].add(name)

        stats["total"] += 1

        # Track morphotype distribution
        morph = row["morphotype"] or "unknown"
        stats["morphotypes"][morph] = stats["morphotypes"].get(morph, 0) + 1

        # Download EM image
        em_url = row["em_file"]
        if em_url:
            stats["has_em"] += 1
            if not args.skip_images:
                ext = Path(em_url).suffix or ".jpg"
                em_path = em_dir / f"{name}{ext}"
                if em_path.exists():
                    row["em_local_path"] = str(em_path)
                    stats["em_downloaded"] += 1
                elif download_image(em_url, em_path):
                    row["em_local_path"] = str(em_path)
                    stats["em_downloaded"] += 1
                    print(f"    EM: {name}{ext}")
                else:
                    print(f"    [FAIL] EM: {name}")

        # Download plaque image
        plaque_url = row["plaque_file"]
        if plaque_url:
            stats["has_plaque"] += 1
            if not args.skip_images:
                ext = Path(plaque_url).suffix or ".jpg"
                plaque_path = plaque_dir / f"{name}{ext}"
                if plaque_path.exists():
                    row["plaque_local_path"] = str(plaque_path)
                    stats["plaque_downloaded"] += 1
                elif download_image(plaque_url, plaque_path):
                    row["plaque_local_path"] = str(plaque_path)
                    stats["plaque_downloaded"] += 1
                    print(f"    Plaque: {name}{ext}")
                else:
                    print(f"    [FAIL] Plaque: {name}")

        # Download EM thumbnail
        em_thumb_url = row["em_thumb_file"]
        if em_thumb_url:
            stats["has_em_thumb"] += 1
            if not args.skip_images:
                ext = Path(em_thumb_url).suffix or ".jpg"
                em_thumb_path = em_thumb_dir / f"{name}{ext}"
                if em_thumb_path.exists():
                    row["em_thumb_local_path"] = str(em_thumb_path)
                    stats["em_thumb_downloaded"] += 1
                elif download_image(em_thumb_url, em_thumb_path):
                    row["em_thumb_local_path"] = str(em_thumb_path)
                    stats["em_thumb_downloaded"] += 1
                    print(f"    EM Thumb: {name}{ext}")
                else:
                    print(f"    [FAIL] EM Thumb: {name}")

        # Download plaque thumbnail
        plaque_thumb_url = row["plaque_thumb_file"]
        if plaque_thumb_url:
            stats["has_plaque_thumb"] += 1
            if not args.skip_images:
                ext = Path(plaque_thumb_url).suffix or ".jpg"
                plaque_thumb_path = plaque_thumb_dir / f"{name}{ext}"
                if plaque_thumb_path.exists():
                    row["plaque_thumb_local_path"] = str(plaque_thumb_path)
                    stats["plaque_thumb_downloaded"] += 1
                elif download_image(plaque_thumb_url, plaque_thumb_path):
                    row["plaque_thumb_local_path"] = str(plaque_thumb_path)
                    stats["plaque_thumb_downloaded"] += 1
                    print(f"    Plaque Thumb: {name}{ext}")
                else:
                    print(f"    [FAIL] Plaque Thumb: {name}")

        # Download FASTA genome sequence
        fasta_url = row["fasta_file"]
        if fasta_url:
            stats["has_fasta"] += 1
            if not args.skip_images:
                ext = Path(fasta_url).suffix or ".fasta"
                fasta_path = fasta_dir / f"{name}{ext}"
                if fasta_path.exists():
                    row["fasta_local_path"] = str(fasta_path)
                    stats["fasta_downloaded"] += 1
                elif download_image(fasta_url, fasta_path):
                    row["fasta_local_path"] = str(fasta_path)
                    stats["fasta_downloaded"] += 1
                    print(f"    FASTA: {name}{ext}")
                else:
                    print(f"    [FAIL] FASTA: {name}")

        # Write row incrementally
        csv_writer.writerow(row)
        csv_file.flush()

        # Progress update every 50 phages
        if stats["total"] % 50 == 0:
            print(f"\n--- Progress: {stats['total']} phages processed ---\n")
            # Write interim summary so progress is saved
            _write_summary(out_dir / "summary.json", stats, args.cluster)

    csv_file.close()

    # Write final summary
    summary_path = out_dir / "summary.json"
    _write_summary(summary_path, stats, args.cluster)

    # Print summary
    print("\n" + "=" * 60)
    print("DOWNLOAD COMPLETE")
    print("=" * 60)
    print(f"Total phages:         {stats['total']}")
    print(f"With EM micrograph:   {stats['has_em']}")
    print(f"With plaque image:    {stats['has_plaque']}")
    print(f"With EM thumbnail:    {stats['has_em_thumb']}")
    print(f"With plaque thumbnail:{stats['has_plaque_thumb']}")
    print(f"With FASTA:           {stats['has_fasta']}")
    if not args.skip_images:
        print(f"EM downloaded:        {stats['em_downloaded']}")
        print(f"Plaque downloaded:    {stats['plaque_downloaded']}")
        print(f"EM thumb downloaded:  {stats['em_thumb_downloaded']}")
        print(f"Plaque thumb downloaded:{stats['plaque_thumb_downloaded']}")
        print(f"FASTA downloaded:     {stats['fasta_downloaded']}")
    print(f"\nMorphotype distribution:")
    for morph, count in sorted(stats["morphotypes"].items()):
        print(f"  {morph:10s}: {count}")
    print(f"\nMetadata CSV: {csv_path}")
    print(f"Summary JSON: {summary_path}")


if __name__ == "__main__":
    main()
