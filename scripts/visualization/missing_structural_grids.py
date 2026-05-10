"""
EM grids for phages missing a structural PHROG class.

Generates one figure per missing class, showing every phage in that group
that has an EM image (clean thumbnail preferred, raw fallback). Used as a
sanity-check that the morphology in these annotation-incomplete genomes
still looks like a phage and is worth feeding to FACE.

Outputs:
    figures/missing_head_em.png
    figures/missing_tail_em.png
"""

import math
from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image

THUMB_DIR = Path("data/processed_data/em_thumb_clean")
RAW_DIR = Path("data/downloaded_data/phagesdb_dataset/em_images")
OUT_DIR = Path("figures")

GROUPS = {
    "missing_head": (
        "Phages missing 'head and packaging' annotation",
        ["Araxxi", "Moyashi"],
    ),
    "missing_tail": (
        "Phages missing 'tail' annotation",
        ["Biozilla", "CrunchyBoi", "Dennebes", "Fabian", "FlowerPower",
         "GalaxyEater", "Geostin", "Gremlin23", "Noely", "Pabst",
         "PigHeaded", "PineapplePluto", "RetrieverFever", "TicTac",
         "Vorvolakos"],
    ),
}


def find_image(name):
    thumb = THUMB_DIR / f"{name}.png"
    if thumb.is_file():
        return thumb, "thumb"
    for p in sorted(RAW_DIR.glob(f"{name}.*")):
        return p, "raw"
    return None, None


def plot_grid(title, names, out_path):
    n = len(names)
    cols = min(5, n)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.4, rows * 2.6))
    axes = [axes] if n == 1 else (axes.flatten() if rows * cols > 1 else [axes])

    for ax, name in zip(axes, names):
        path, kind = find_image(name)
        if path is None:
            ax.text(0.5, 0.5, "no image", ha="center", va="center")
            ax.set_facecolor("#eee")
        else:
            img = Image.open(path).convert("RGB")
            ax.imshow(img)
        label = name if kind != "raw" else f"{name} (raw)"
        ax.set_title(label, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}  ({n} phages)")


def main():
    for key, (title, names) in GROUPS.items():
        plot_grid(title, names, OUT_DIR / f"{key}_em.png")


if __name__ == "__main__":
    main()
