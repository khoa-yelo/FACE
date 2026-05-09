"""Standardize PhagesDB thumbnail images for ML.

Input:
  data/downloaded_data/phagesdb_dataset/em_thumb_images/        (mixed format/mode/size)
  data/downloaded_data/phagesdb_dataset/plaque_thumb_images/    (mixed format/mode/size)

Output:
  data/processed_data/em_thumb_clean/<name>.png        200x200 RGB, lossless
  data/processed_data/plaque_thumb_clean/<name>.png    200x200 RGB, lossless
  data/processed_data/standardize_failed.tsv           (input, error)

Per-image transformations:
  1. Multi-frame (MPO/GIF): take frame 0
  2. 16-bit modes (I;16, I;16B): scale to 8-bit
  3. Alpha (RGBA/LA/P-with-transparency): composite onto white
  4. Convert mode to RGB
  5. Resize to 200x200 with Lanczos if needed
  6. Save as PNG (lossless)

Idempotent: skips outputs that already exist with nonzero size.
"""
import argparse
import csv
import multiprocessing as mp
import sys
from pathlib import Path

from PIL import Image, ImageFile

# Allow large EM images that PIL flags as decompression bombs (e.g. CooJo.jpg)
Image.MAX_IMAGE_PIXELS = None
# Be tolerant of truncated/partial files
ImageFile.LOAD_TRUNCATED_IMAGES = True

INPUT_ROOT = Path("data/downloaded_data/phagesdb_dataset")
OUTPUT_ROOT = Path("data/processed_data")
TARGET_SIZE = (200, 200)


def standardize_one(args):
    """Worker: returns (input_path, status, message)."""
    src, dst = args
    try:
        if dst.exists() and dst.stat().st_size > 0:
            return (str(src), "skip", "already exists")

        with Image.open(src) as im:
            # 1) multi-frame: take first frame
            if getattr(im, "n_frames", 1) > 1:
                im.seek(0)

            # 2) 16-bit -> 8-bit (numpy path; PIL.point() doesn't handle I;16B)
            if im.mode in ("I;16", "I;16B", "I;16L", "I"):
                import numpy as np
                arr = np.asarray(im, dtype=np.uint16 if im.mode != "I" else np.int32)
                arr8 = (arr.astype(np.float32) / 256.0).clip(0, 255).astype(np.uint8)
                im = Image.fromarray(arr8, mode="L")

            # 3) alpha compositing onto white
            if im.mode in ("RGBA", "LA"):
                bg = Image.new("RGB", im.size, (255, 255, 255))
                # split alpha as the mask
                alpha = im.getchannel("A")
                im = im.convert("RGBA")
                bg.paste(im, mask=alpha)
                im = bg
            elif im.mode == "P" and "transparency" in im.info:
                im = im.convert("RGBA")
                bg = Image.new("RGB", im.size, (255, 255, 255))
                alpha = im.getchannel("A")
                bg.paste(im, mask=alpha)
                im = bg

            # 4) ensure RGB
            if im.mode != "RGB":
                im = im.convert("RGB")

            # 5) resize if needed
            if im.size != TARGET_SIZE:
                im = im.resize(TARGET_SIZE, Image.LANCZOS)

            # 6) save
            dst.parent.mkdir(parents=True, exist_ok=True)
            im.save(dst, format="PNG", optimize=True)
        return (str(src), "ok", "")
    except Exception as e:
        return (str(src), "fail", repr(e))


def collect_jobs(in_dir: Path, out_dir: Path):
    jobs = []
    for src in sorted(in_dir.iterdir()):
        if not src.is_file():
            continue
        dst = out_dir / (src.stem + ".png")
        jobs.append((src, dst))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None,
                    help="Process only first N images (for testing)")
    args = ap.parse_args()

    pairs = [
        ("em_thumb_images", "em_thumb_clean"),
        ("plaque_thumb_images", "plaque_thumb_clean"),
    ]

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    failed_path = OUTPUT_ROOT / "standardize_failed.tsv"
    with open(failed_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["input", "error"])

        for sub_in, sub_out in pairs:
            in_dir = INPUT_ROOT / sub_in
            out_dir = OUTPUT_ROOT / sub_out
            jobs = collect_jobs(in_dir, out_dir)
            if args.limit:
                jobs = jobs[: args.limit]

            print(f"\n=== {sub_in} -> {sub_out}: {len(jobs)} images ===")
            ok = skipped = failed = 0

            with mp.Pool(processes=args.workers) as pool:
                for i, (src, status, msg) in enumerate(
                    pool.imap_unordered(standardize_one, jobs, chunksize=32),
                    start=1,
                ):
                    if status == "ok":
                        ok += 1
                    elif status == "skip":
                        skipped += 1
                    else:
                        failed += 1
                        writer.writerow([src, msg])
                    if i % 1000 == 0:
                        print(f"  {i}/{len(jobs)}  ok={ok} skip={skipped} fail={failed}")

            print(f"  done: ok={ok} skip={skipped} fail={failed}")

    print(f"\nFailures logged to: {failed_path}")


if __name__ == "__main__":
    main()
