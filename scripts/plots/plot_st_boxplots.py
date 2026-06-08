"""Per-ST gene-load boxplots colored by dominant phylogroup.

For each Sequence Type with >= MIN_ST_SIZE isolates, plot the distribution of
focal-gene counts per isolate. STs are colored by their dominant phylogroup.
A right-hand sidebar reports the per-ST isolate count on a log scale.
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

from utils import save_fig


MIN_ST_SIZE = 10  # min isolates per ST to include

COUNT_COL = "n_genes"
COUNT_LABEL = "tot n. of genes in isolate"

PG_PALETTE = {
    "A": "#66c2a5",
    "B1": "#fc8d62",
    "B2": "#8da0cb",
    "C": "#b3b3b3",
    "D": "#e78ac3",
    "E": "#a6d854",
    "F": "#ffd92f",
    "Shigella": "#999999",
}
DEFAULT_PG_COLOR = "lightgray"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pa", required=True, help="Focal-gene presence/absence CSV")
    p.add_argument("--metadata", required=True, help="Horesh F1 genome metadata CSV")
    p.add_argument("--out-boxplot", required=True, nargs="+", help="Output path(s) for the boxplot+counts figure (format inferred from extension; pass one per format)")
    p.add_argument("--out-csv", required=True, help="Output per-isolate gene-count CSV")
    p.add_argument("--high-load-threshold", type=int, required=True, help="Gene-count threshold marking 'high-load' isolates")
    return p.parse_args()


def load_data(args):
    """Read and preprocess inputs; return (pa, meta)."""
    pa = pd.read_csv(args.pa, index_col=0)
    pa.index = pa.index.astype(int)

    meta = pd.read_csv(args.metadata)
    meta = meta.dropna(subset=["name_in_presence_absence"])
    meta["Phylogroup"] = meta["Phylogroup"].replace({"Not determined": "Not Determined"})

    return pa, meta


def compute_isolate_counts(pa, meta):
    """Per-isolate focal-gene count with ST and phylogroup, over all pa columns."""
    strain_meta = meta.set_index("name_in_presence_absence")
    counts = pd.DataFrame({COUNT_COL: pa.sum(axis=0)})
    counts["ST"] = strain_meta["ST"].reindex(counts.index)
    counts["Phylogroup"] = strain_meta["Phylogroup"].reindex(counts.index)
    counts.index.name = "isolate"
    print(f"Isolates: {len(counts)}")
    return counts


def prepare_plot_data(isolate_counts, meta):
    """Filter to valid STs and build plot ordering/palette.

    Returns (count_df, st_order, palette, st_to_pg).
    """
    count_df = isolate_counts.dropna(subset=["ST"]).copy()
    count_df = count_df[~count_df["ST"].astype(str).str.endswith("~")]

    st_counts = count_df["ST"].value_counts()
    valid_sts = st_counts[st_counts >= MIN_ST_SIZE].index
    count_df = count_df[count_df["ST"].isin(valid_sts)]

    print(f"Strains after filtering: {len(count_df)}")
    print(f"Valid STs (>={MIN_ST_SIZE} members): {len(valid_sts)}")

    st_to_pg = {}
    for st in valid_sts:
        pg = meta[meta["ST"].astype(str) == str(st)]["Phylogroup"].mode()
        st_to_pg[st] = pg.iloc[0] if len(pg) > 0 else "Not Determined"

    st_order = count_df.groupby("ST")[COUNT_COL].median().sort_values().index
    palette = {st: PG_PALETTE.get(st_to_pg.get(st, "?"), DEFAULT_PG_COLOR) for st in st_order}

    return count_df, st_order, palette, st_to_pg


def draw_st_boxplot(ax, count_df, col, label, st_order, palette, st_to_pg, high_load_threshold):
    sns.boxplot(
        y="ST",
        x=col,
        hue="ST",
        data=count_df,
        order=st_order,
        palette=palette,
        showfliers=False,
        legend=False,
        ax=ax,
    )
    n_collections_before = len(ax.collections)
    sns.stripplot(
        y="ST",
        x=col,
        data=count_df,
        order=st_order,
        color="k",
        marker=".",
        alpha=0.3,
        size=2,
        ax=ax,
    )
    for c in ax.collections[n_collections_before:]:
        c.set_rasterized(True)
    ax.set_xlabel(label, fontsize=12)
    ax.set_ylabel("Sequence type (ST)", fontsize=12)
    ax.axvline(high_load_threshold, color="k", linestyle="--")
    ax.text(
        high_load_threshold + 5,
        1,
        f"high-load threshold ({high_load_threshold} genes)",
        va="top",
        ha="left",
        rotation=90,
        fontsize=12,
    )
    present_pgs = {st_to_pg[st] for st in st_order}
    legend_handles = [
        Patch(color=PG_PALETTE.get(pg, DEFAULT_PG_COLOR), label=pg)
        for pg in list(PG_PALETTE.keys()) + ["Not Determined"]
        if pg in present_pgs
    ]
    ax.legend(
        handles=legend_handles,
        title="Phylogroup",
        loc="lower left",
        frameon=True,
        fontsize=11,
        title_fontsize=12,
    )


def plot_st_boxplots(count_df, st_order, palette, st_to_pg, high_load_threshold, out_boxplot):
    """Render the per-ST boxplot + isolate-count sidebar and save to `out_boxplot`."""
    fig, (ax_box, ax_bar) = plt.subplots(
        1,
        2,
        figsize=(8, 10),
        sharey=True,
        gridspec_kw={"width_ratios": [4, 1], "wspace": 0.05},
    )
    draw_st_boxplot(ax_box, count_df, COUNT_COL, COUNT_LABEL, st_order, palette, st_to_pg, high_load_threshold)
    st_sizes = count_df["ST"].value_counts().reindex(st_order)
    ax_bar.barh(
        range(len(st_order)),
        st_sizes.values,
        color="lightgray",
        edgecolor="gray",
        linewidth=0.4,
    )
    ax_bar.set_xscale("log")
    ax_bar.set_xlabel("# isolates", fontsize=12)
    ax_bar.tick_params(labelleft=False)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.set_ylim(len(st_order) - 0.5, -0.5)
    plt.tight_layout()
    save_fig(fig, out_boxplot)
    print(f"Saved {out_boxplot}")


def main():
    args = parse_args()

    pa, meta = load_data(args)

    isolate_counts = compute_isolate_counts(pa, meta)
    isolate_counts.to_csv(args.out_csv)
    print(f"Saved {args.out_csv}")

    count_df, st_order, palette, st_to_pg = prepare_plot_data(isolate_counts, meta)
    plot_st_boxplots(count_df, st_order, palette, st_to_pg, args.high_load_threshold, args.out_boxplot)


if __name__ == "__main__":
    main()
