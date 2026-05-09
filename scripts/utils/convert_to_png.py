#!/usr/bin/env python3
"""Convert all images in an input folder to PNG format in an output folder."""

import argparse
import os
from PIL import Image


def convert_images_to_png(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    supported_extensions = {".jpg", ".jpeg", ".bmp", ".gif", ".tiff", ".tif", ".webp", ".ico"}
    converted = 0
    skipped = 0

    for filename in os.listdir(input_folder):
        filepath = os.path.join(input_folder, filename)
        if not os.path.isfile(filepath):
            continue

        name, ext = os.path.splitext(filename)
        ext_lower = ext.lower()

        if ext_lower == ".png":
            # Already PNG, just copy
            img = Image.open(filepath)
            img.save(os.path.join(output_folder, filename))
            converted += 1
            print(f"Copied: {filename}")
            continue

        if ext_lower not in supported_extensions:
            print(f"Skipped (unsupported): {filename}")
            skipped += 1
            continue

        try:
            img = Image.open(filepath)
            img = img.convert("RGBA") if img.mode in ("RGBA", "LA", "PA") else img.convert("RGB")
            out_path = os.path.join(output_folder, f"{name}.png")
            img.save(out_path, "PNG")
            converted += 1
            print(f"Converted: {filename} -> {name}.png")
        except Exception as e:
            print(f"Error converting {filename}: {e}")
            skipped += 1

    print(f"\nDone. Converted: {converted}, Skipped: {skipped}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert all images in a folder to PNG.")
    parser.add_argument("input_folder", help="Path to input folder containing images")
    parser.add_argument("output_folder", help="Path to output folder for PNG images")
    args = parser.parse_args()

    convert_images_to_png(args.input_folder, args.output_folder)
