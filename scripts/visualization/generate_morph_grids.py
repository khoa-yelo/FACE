import csv
import os
import random
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METADATA_CSV = 'data/processed_data/metadata.csv'
OUT_DIR = 'figures/morph_examples'
os.makedirs(OUT_DIR, exist_ok=True)

# Collect images per morphotype
morphs = {}
with open(METADATA_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        m = row['morphotype'].strip()
        if not m:
            m = 'Unknown'
        thumb = row.get('em_thumb_local_path', '').strip()
        if thumb:
            full_path = thumb
            if os.path.exists(full_path):
                morphs.setdefault(m, []).append((row['phage_name'], full_path))

# For each morphotype, pick up to 9 random examples and make a 3x3 grid
for morph_name, items in morphs.items():
    n_samples = min(9, len(items))
    samples = random.sample(items, n_samples)

    ncols = 3
    nrows = (n_samples + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows))
    fig.suptitle(f'Morphotype: {morph_name} ({len(items)} total images)',
                 fontsize=18, fontweight='bold', y=1.02)

    if nrows == 1:
        axes = [axes]

    for i in range(nrows * ncols):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        if i < n_samples:
            name, path = samples[i]
            try:
                img = Image.open(path).convert('RGB')
                ax.imshow(img)
                ax.set_title(name, fontsize=11)
            except Exception as e:
                ax.text(0.5, 0.5, f'Error loading\n{name}', ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')

    plt.tight_layout()
    safe_name = morph_name.replace(' ', '_')
    out_path = os.path.join(OUT_DIR, f'morph_{safe_name}.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_path}')

print('\nDone! All grid images saved to', OUT_DIR)
