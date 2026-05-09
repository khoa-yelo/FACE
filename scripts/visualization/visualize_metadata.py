#!/usr/bin/env python3
"""Visualize PhagesDB metadata breakdown by morphology, cluster, sequencing status, etc."""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


def load_data(csv_path):
    df = pd.read_csv(csv_path)
    df["morphotype"] = df["morphotype"].fillna("Unknown").replace("", "Unknown")
    df["cluster"] = df["cluster"].fillna("Singleton").replace("", "Singleton")
    df["temperate"] = df["temperate"].map({True: "Temperate", False: "Lytic", "True": "Temperate", "False": "Lytic"}).fillna("Unknown")
    df["seq_status"] = df["seq_finished"].map({True: "Sequenced", False: "Not Sequenced", "True": "Sequenced", "False": "Not Sequenced"}).fillna("Unknown")
    df["has_em"] = df["em_file"].notna() & (df["em_file"] != "")
    df["has_plaque"] = df["plaque_file"].notna() & (df["plaque_file"] != "")
    return df


def plot_morphotype_em_thumb(df, ax):
    df = df.copy()
    df["has_em_thumb"] = df["em_thumb_file"].notna() & (df["em_thumb_file"] != "")
    ct = pd.crosstab(df["morphotype"], df["has_em_thumb"])
    ct.columns = ["No EM Thumb", "Has EM Thumb"]
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    ct.plot.bar(ax=ax, stacked=True, color=["#e74c3c", "#2ecc71"], edgecolor="black")
    ax.set_title("Morphotype: EM Thumbnail Availability")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title="")
    for i, (idx, row) in enumerate(ct.iterrows()):
        total = row.sum()
        ax.text(i, total + 100, f"{total:,}", ha='center', fontsize=7)


def plot_top_clusters(df, ax, top_n=20):
    non_singleton = df[df["cluster"] != "Singleton"]
    ct = pd.crosstab(non_singleton["cluster"], non_singleton["morphotype"])
    top_clusters = non_singleton["cluster"].value_counts().head(top_n).index
    ct = ct.loc[top_clusters]
    morph_order = ["SIPHO", "MYO", "PODO", "TECTI", "Unknown"]
    morph_order = [m for m in morph_order if m in ct.columns]
    ct = ct[morph_order]
    color_map = {"SIPHO": "#3498db", "MYO": "#e74c3c", "PODO": "#2ecc71", "TECTI": "#9b59b6", "Unknown": "#95a5a6"}
    colors = [color_map.get(m, "#95a5a6") for m in morph_order]
    ct.plot.bar(ax=ax, stacked=True, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_title(f"Top {top_n} Clusters by Morphotype (excl. Singletons)")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title="", fontsize=7)
    for i, (idx, row) in enumerate(ct.iterrows()):
        total = row.sum()
        ax.text(i, total + 2, f"{total:,}", ha='center', fontsize=6)


def plot_singleton(df, ax):
    singleton = df[df["cluster"] == "Singleton"]
    ct = singleton["morphotype"].value_counts()
    color_map = {"SIPHO": "#3498db", "MYO": "#e74c3c", "PODO": "#2ecc71", "TECTI": "#9b59b6", "Unknown": "#95a5a6"}
    colors = [color_map.get(m, "#95a5a6") for m in ct.index]
    bars = ax.bar(ct.index, ct.values, color=colors, edgecolor="black")
    ax.set_title(f"Singletons by Morphotype (n={len(singleton):,})")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=45)
    for bar, v in zip(bars, ct.values):
        ax.text(bar.get_x() + bar.get_width()/2, v + 20, f"{v:,}", ha='center', fontsize=7)


def plot_temperate(df, ax):
    ct = pd.crosstab(df["morphotype"], df["temperate"])
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    col_order = [c for c in ["Temperate", "Lytic", "Unknown"] if c in ct.columns]
    ct = ct[col_order]
    color_map = {"Temperate": "#3498db", "Lytic": "#e67e22", "Unknown": "#95a5a6"}
    colors = [color_map[c] for c in col_order]
    ct.plot.bar(ax=ax, stacked=False, color=colors, edgecolor="black")
    ax.set_title("Temperate vs Lytic by Morphotype")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title="")
    for container in ax.containers:
        ax.bar_label(container, fontsize=7, padding=2)


def plot_morphotype_seq(df, ax):
    ct = pd.crosstab(df["morphotype"], df["seq_status"])
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    col_order = [c for c in ["Sequenced", "Not Sequenced"] if c in ct.columns]
    ct = ct[col_order]
    colors = ["#2ecc71", "#e74c3c"]
    ct.plot.bar(ax=ax, stacked=True, color=colors, edgecolor="black")
    ax.set_title("Morphotype by Sequencing Status")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title="")
    for i, (idx, row) in enumerate(ct.iterrows()):
        total = row.sum()
        ax.text(i, total + 100, f"{total:,}", ha='center', fontsize=7)


def plot_paired_seq_em(df, ax):
    df = df.copy()
    df["has_em"] = df["em_file"].notna() & (df["em_file"] != "")
    df["has_seq"] = df["seq_finished"].isin([True, "True"])
    df["paired"] = df["has_em"] & df["has_seq"]

    morpho_order = df["morphotype"].value_counts().index
    result = pd.DataFrame({
        "Has Sequence": df.groupby("morphotype")["has_seq"].sum(),
        "Has EM Image": df.groupby("morphotype")["has_em"].sum(),
        "Paired (Seq + EM)": df.groupby("morphotype")["paired"].sum(),
    }).loc[morpho_order]

    result.plot.bar(ax=ax, edgecolor="black")
    ax.set_title("Paired Sequence + EM Image by Morphotype")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title="", fontsize=8)
    for container in ax.containers:
        ax.bar_label(container, fontsize=7, padding=2)


def main():
    parser = argparse.ArgumentParser(description="Visualize PhagesDB metadata")
    parser.add_argument("--csv", default="data/processed_data/metadata.csv",
                        help="Path to metadata.csv")
    parser.add_argument("--output", default="figures/phagesdb_metadata_summary.png",
                        help="Output figure path")
    args = parser.parse_args()

    df = load_data(args.csv)

    print(f"Total phages: {len(df)}")
    print(f"Morphotypes: {df['morphotype'].value_counts().to_dict()}")
    print(f"Sequencing: {df['seq_status'].value_counts().to_dict()}")
    print(f"Temperate: {df['temperate'].value_counts().to_dict()}")
    print(f"Clusters: {df['cluster'].nunique()}")
    print(f"With EM: {df['has_em'].sum()}, With Plaque: {df['has_plaque'].sum()}")

    fig, axes = plt.subplots(3, 3, figsize=(20, 16))
    fig.suptitle(f"PhagesDB Dataset Overview (n={len(df):,})", fontsize=14, fontweight='bold')

    plot_morphotype_em_thumb(df, axes[0, 0])
    plot_top_clusters(df, axes[0, 1])
    plot_singleton(df, axes[0, 2])
    plot_temperate(df, axes[1, 0])
    plot_morphotype_seq(df, axes[1, 1])
    plot_paired_seq_em(df, axes[1, 2])
    axes[2, 0].axis("off")
    axes[2, 1].axis("off")
    axes[2, 2].axis("off")

    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to: {args.output}")


if __name__ == "__main__":
    main()
