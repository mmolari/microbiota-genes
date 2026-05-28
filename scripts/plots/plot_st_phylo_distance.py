"""ST-level phylogenetic distances from the Horesh tree.

Computes pairwise leaf distances on tree_500, aggregates them to ST-level
mean pairwise distances, and renders three figures:
  - per-ST isolate counts (total vs in-tree)
  - clustered heatmap of ST-level mean pairwise distances
  - the tree itself with phylogroup-coloured tips and ST labels
Also writes two derived CSVs (isolate counts and pairwise-distance matrix).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from utils import save_fig


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

MIN_ISOLATES = 10  # minimum isolates per ST to include in distance matrix


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree", required=True, help="Horesh tree_500 Newick")
    p.add_argument("--metadata", required=True, help="Horesh F1 genome metadata CSV")
    p.add_argument(
        "--out-isolates",
        required=True,
        nargs="+",
        help="Output path(s) for the isolates-per-ST barplot (format inferred from extension)",
    )
    p.add_argument(
        "--out-distance",
        required=True,
        nargs="+",
        help="Output path(s) for the ST pairwise-distance heatmap",
    )
    p.add_argument(
        "--out-counts-csv", required=True, help="Output CSV with per-ST isolate counts"
    )
    p.add_argument(
        "--out-distance-csv",
        required=True,
        help="Output CSV with ST pairwise mean-distance matrix",
    )
    return p.parse_args()


def root_paths(tree):
    paths = {}
    for leaf in tree.get_terminals():
        cum = 0.0
        d = {id(tree.root): 0.0}
        for c in tree.get_path(leaf):
            cum += c.branch_length or 0.0
            d[id(c)] = cum
        paths[leaf.name] = d
    return paths


def main():
    args = parse_args()

    tree = Phylo.read(args.tree, "newick")

    meta = pd.read_csv(args.metadata)
    meta = meta.dropna(subset=["name_in_presence_absence"]).set_index(
        "name_in_presence_absence"
    )
    meta["ST"] = meta["ST"].astype(str).str.strip()
    meta["Phylogroup"] = meta["Phylogroup"].replace(
        {"Not determined": "Not Determined"}
    )
    meta = meta[~meta["ST"].str.endswith("~")]

    leaves_all = [t.name for t in tree.get_terminals()]
    leaves = [l for l in leaves_all if l in meta.index]
    leaf_st = meta.loc[leaves, "ST"]
    print(f"tree leaves: {len(leaves_all)}, kept after ~ filter: {len(leaves)}")

    # Per-ST counts
    st_in_tree = leaf_st.value_counts()
    st_total = meta["ST"].value_counts()
    counts = (
        pd.DataFrame({"n_in_tree": st_in_tree, "n_total": st_total})
        .fillna(0)
        .astype(int)
    )
    counts["fraction_in_tree"] = counts["n_in_tree"] / counts["n_total"]
    counts = counts.sort_values("n_total", ascending=False)
    counts.index.name = "ST"
    counts.to_csv(args.out_counts_csv)

    # Barplot: isolates per ST
    plot_df = counts[counts["n_total"] >= 3].copy()
    fig, ax = plt.subplots(figsize=(max(8, 0.13 * len(plot_df)), 5))
    x = np.arange(len(plot_df))
    ax.bar(x, plot_df["n_total"], color="lightgray", label="total in dataset")
    ax.bar(x, plot_df["n_in_tree"], color="steelblue", label="in tree_500")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df.index, rotation=90, fontsize=6)
    ax.set_xlabel("ST")
    ax.axhline(
        MIN_ISOLATES,
        color="k",
        linestyle="--",
        label=f"{MIN_ISOLATES} isolates",
    )
    ax.set_ylabel("number of isolates")
    ax.legend()
    ax.set_xlim(-0.5, len(plot_df) - 0.5)
    fig.tight_layout()
    save_fig(fig, args.out_isolates)

    # Pairwise leaf distances via shared root paths
    paths = root_paths(tree)
    leaf_clade = {t.name: t for t in tree.get_terminals()}
    N = len(leaves)
    D = np.zeros((N, N))
    for i, a in enumerate(leaves):
        pa_a = paths[a]
        da = pa_a[id(leaf_clade[a])]
        for j in range(i + 1, N):
            b = leaves[j]
            pb = paths[b]
            common = pa_a.keys() & pb.keys()
            lca_depth = max(pa_a[c] for c in common)
            D[i, j] = D[j, i] = (da - lca_depth) + (pb[id(leaf_clade[b])] - lca_depth)

    # Aggregate to ST-level mean pairwise distances
    keep_sts = counts[
        (counts["n_total"] >= MIN_ISOLATES) & (counts["n_in_tree"] > 0)
    ].index.tolist()
    leaf_to_st = leaf_st.to_dict()
    st_idx = {
        s: [i for i, l in enumerate(leaves) if leaf_to_st[l] == s] for s in keep_sts
    }

    M = pd.DataFrame(index=keep_sts, columns=keep_sts, dtype=float)
    for s1 in keep_sts:
        for s2 in keep_sts:
            block = D[np.ix_(st_idx[s1], st_idx[s2])]
            if s1 == s2:
                iu = np.triu_indices_from(block, k=1)
                M.loc[s1, s2] = block[iu].mean() if iu[0].size else 0.0
            else:
                M.loc[s1, s2] = block.mean()

    # Hierarchical clustering order
    condensed = squareform(M.values, checks=False)
    order_idx = leaves_list(linkage(condensed, method="average"))
    M_ord = M.iloc[order_idx, order_idx]
    M_ord.to_csv(args.out_distance_csv)

    # ST -> dominant phylogroup
    st_to_pg = {}
    for st in keep_sts:
        pg = meta.loc[meta["ST"] == st, "Phylogroup"].mode()
        st_to_pg[st] = pg.iloc[0] if len(pg) > 0 else "Not Determined"

    # Distance heatmap
    fig, ax = plt.subplots(
        figsize=(max(5, 0.18 * len(M_ord)), max(4, 0.18 * len(M_ord)))
    )
    im = ax.imshow(M_ord.values, cmap="viridis", aspect="equal")
    ax.set_xticks(range(len(M_ord)))
    ax.set_xticklabels(M_ord.index, rotation=90, fontsize=10)
    ax.set_yticks(range(len(M_ord)))
    ax.set_yticklabels(M_ord.index, fontsize=10)
    ax.set_xlabel("Sequence Type")
    ax.set_ylabel("Sequence Type")
    for tick, st in zip(ax.get_xticklabels(), M_ord.index):
        tick.set_color(PG_PALETTE.get(st_to_pg.get(st), DEFAULT_PG_COLOR))
        tick.set_fontweight("bold")
    for tick, st in zip(ax.get_yticklabels(), M_ord.index):
        tick.set_color(PG_PALETTE.get(st_to_pg.get(st), DEFAULT_PG_COLOR))
        tick.set_fontweight("bold")
    fig.colorbar(
        im, ax=ax, label="mean pairwise distance", shrink=0.46, fraction=0.03, pad=0.02
    )

    present_pgs = {st_to_pg[st] for st in M_ord.index}
    legend_handles = [
        Patch(color=PG_PALETTE.get(pg, DEFAULT_PG_COLOR), label=pg)
        for pg in list(PG_PALETTE.keys()) + ["Not Determined"]
        if pg in present_pgs
    ]
    ax.legend(
        handles=legend_handles,
        title="Phylogroup",
        loc="upper left",
        bbox_to_anchor=(1.15, 1.0),
        frameon=True,
        fontsize=9.2,
        title_fontsize=10.35,
    )
    fig.tight_layout()
    save_fig(fig, args.out_distance)


if __name__ == "__main__":
    main()
