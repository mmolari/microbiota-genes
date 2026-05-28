"""Phylogenetic tree + focal-gene presence/absence heatmap.

Combines:
  - the Horesh tree_500 pruned to strains in the focal-gene PA matrix
  - a strip annotating each tip's ST (top STs by isolate count are coloured)
  - a clustered heatmap of focal-gene PA (non-core genes only)
  - a marginal frequency bar (gene frequency in the full dataset)
  - a highlight rectangle around a specified ST x gene-id range block.
Also exports a supplementary CSV of the rendered PA sub-matrix.
"""

import argparse
import copy
import csv
import pathlib
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle
from scipy.cluster.hierarchy import linkage, to_tree

from utils import save_fig


CORE_FRACTION = 0.95     # genes with freq >= this are treated as core and excluded
ST_MIN_STRAINS = 5       # min strains per ST to colour-highlight on the ST strip
HIGHLIGHT_ST = "131"     # ST used to draw the highlight rectangle
HIGHLIGHT_GENE_A = 3258  # yehH — one corner of the highlight rectangle
HIGHLIGHT_GENE_B = 4518  # other corner of the highlight rectangle


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pa", required=True, help="Focal-gene presence/absence CSV")
    p.add_argument("--gene-info", required=True, help="Focal gene_info CSV")
    p.add_argument("--metadata", required=True, help="Horesh F1 genome metadata CSV")
    p.add_argument("--tree", required=True, help="Horesh tree_500 Newick")
    p.add_argument("--out-heatmap", required=True, nargs="+", help="Output path(s) for the heatmap figure (format inferred from extension)")
    p.add_argument("--out-pa-csv", required=True, help="Output CSV with the rendered PA sub-matrix")
    return p.parse_args()


def prune_tree(tree, keep_tips):
    tree = copy.deepcopy(tree)
    keep_set = set(keep_tips)
    changed = True
    while changed:
        changed = False
        for clade in list(tree.find_clades(order="postorder")):
            to_remove = [
                c for c in clade.clades if c.is_terminal() and c.name not in keep_set
            ]
            for c in to_remove:
                clade.clades.remove(c)
                changed = True
            to_remove_internal = [
                c for c in clade.clades if not c.is_terminal() and len(c.clades) == 0
            ]
            for c in to_remove_internal:
                clade.clades.remove(c)
                changed = True
            if not clade.is_terminal() and len(clade.clades) == 1:
                child = clade.clades[0]
                bl = (child.branch_length or 0) + (clade.branch_length or 0)
                child.branch_length = bl
                for parent in tree.find_clades():
                    if clade in parent.clades:
                        idx = parent.clades.index(clade)
                        parent.clades[idx] = child
                        changed = True
                        break
    return tree


def tree_segments(tree):
    """Rectangular phylogram with tips at y=0 and root at y=max_depth."""
    terminals = tree.get_terminals()
    depths = tree.depths()
    x_pos = {t: i for i, t in enumerate(terminals)}
    for clade in tree.find_clades(order="postorder"):
        if not clade.is_terminal():
            child_xs = [x_pos[c] for c in clade.clades]
            x_pos[clade] = sum(child_xs) / len(child_xs)
    segs = []
    for clade in tree.find_clades(order="preorder"):
        if not clade.is_terminal():
            child_xs = [x_pos[c] for c in clade.clades]
            segs.append(
                ((min(child_xs), depths[clade]), (max(child_xs), depths[clade]))
            )
        for child in clade.clades:
            segs.append(((x_pos[child], depths[clade]), (x_pos[child], depths[child])))
    return segs, x_pos, depths


def _ordered_leaves(node, freq_arr):
    if node.is_leaf():
        return [node.id]
    left = _ordered_leaves(node.get_left(), freq_arr)
    right = _ordered_leaves(node.get_right(), freq_arr)
    return (
        (left + right)
        if freq_arr[left].mean() >= freq_arr[right].mean()
        else (right + left)
    )


def main():
    args = parse_args()

    pa = pd.read_csv(args.pa, index_col=0)
    pa.index = pa.index.astype(int)

    meta = pd.read_csv(args.metadata)
    meta = meta.dropna(subset=["name_in_presence_absence"])
    meta["ST"] = meta["ST"].astype(str).str.strip()
    strain_st = meta.set_index("name_in_presence_absence")["ST"].to_dict()

    gene_info = pd.read_csv(args.gene_info)
    id_to_name = gene_info.set_index("geneId")["gene_name"].to_dict()

    tree = Phylo.read(args.tree, "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    tree = prune_tree(tree, set(pa.columns))
    tip_order = [t.name for t in tree.get_terminals()]
    tip_order = [t for t in tip_order if t in pa.columns]
    n_strains = len(tip_order)
    print(f"strains in tree after pruning: {n_strains}")

    segs, x_pos, depths = tree_segments(tree)
    max_depth = max(depths.values())

    full_freq = pa.mean(axis=1)
    noncore = full_freq.index[full_freq < CORE_FRACTION]
    print(f"non-core genes (freq < {CORE_FRACTION}): {len(noncore)} / {len(full_freq)}")

    sub_g = pa.loc[noncore, tip_order]
    gene_link = linkage(sub_g.values, method="average", metric="hamming")
    freq_arr = full_freq.loc[noncore].values
    gene_tree = to_tree(gene_link)

    leaves = _ordered_leaves(gene_tree, freq_arr)
    gene_order_clust = sub_g.index[leaves]
    n_genes = len(gene_order_clust)

    # Highlight STs: every ST with >= ST_MIN_STRAINS strains in pruned tree
    st_counts = Counter(strain_st.get(s, "") for s in tip_order)
    highlight_sts = sorted(
        [st for st, n in st_counts.items() if st and n >= ST_MIN_STRAINS],
        key=lambda s: -st_counts[s],
    )
    print(
        f"highlighted STs (>={ST_MIN_STRAINS} strains): "
        f"{len(highlight_sts)} -> {[(st, st_counts[st]) for st in highlight_sts]}"
    )
    qual_colors = (
        list(colormaps["tab20"].colors)
        + list(colormaps["tab20b"].colors)
        + list(colormaps["tab20c"].colors)
    )
    assert len(highlight_sts) <= len(qual_colors), (
        f"too many highlighted STs ({len(highlight_sts)}) for combined palette "
        f"({len(qual_colors)})"
    )
    st_to_color = {st: qual_colors[i] for i, st in enumerate(highlight_sts)}
    st_to_idx = {st: i for i, st in enumerate(highlight_sts)}
    n_other = len(highlight_sts)
    strip_categories = [st_to_idx.get(strain_st.get(s, ""), n_other) for s in tip_order]
    strip_arr = np.array(strip_categories).reshape(1, -1)
    strip_cmap = ListedColormap([st_to_color[st] for st in highlight_sts] + ["white"])

    def gene_label(gid):
        name = id_to_name.get(gid)
        if pd.notna(name) and name not in (None, ""):
            return str(name)
        return f"gene {int(gid)}"

    # Highlight region: HIGHLIGHT_ST columns x (HIGHLIGHT_GENE_A .. HIGHLIGHT_GENE_B) rows
    gene_pos = pd.Index(gene_order_clust)
    y_a = gene_pos.get_loc(HIGHLIGHT_GENE_A)
    y_b = gene_pos.get_loc(HIGHLIGHT_GENE_B)
    y0, y1 = sorted((y_a, y_b))

    hl_st_idx = [i for i, s in enumerate(tip_order) if strain_st.get(s, "") == HIGHLIGHT_ST]
    assert len(hl_st_idx) > 0, f"no ST{HIGHLIGHT_ST} isolates found in tip_order"
    x0, x1 = min(hl_st_idx), max(hl_st_idx)
    print(f"ST{HIGHLIGHT_ST} block: {len(hl_st_idx)} strains in x=[{x0}, {x1}]")
    print(f"gene_a at row {y_a}, gene_b at row {y_b} (highlight rows {y0}..{y1})")

    heat = pa.loc[gene_order_clust, tip_order]

    fig = plt.figure(figsize=(14, 22))
    gs = fig.add_gridspec(
        3,
        2,
        width_ratios=(8, 1),
        height_ratios=(1.2, 0.08, 8),
        hspace=0.02,
        wspace=0.02,
    )
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_legend = fig.add_subplot(gs[0, 1])
    ax_st = fig.add_subplot(gs[1, 0])
    ax_heat = fig.add_subplot(gs[2, 0])
    ax_freq = fig.add_subplot(gs[2, 1])

    ax_heat.imshow(
        heat.values,
        aspect="auto",
        cmap="Greys",
        interpolation="none",
        origin="upper",
        extent=[-0.5, n_strains - 0.5, n_genes, 0],
        vmin=0,
        vmax=1,
    )
    if n_genes > 5:
        ax_heat.hlines(
            np.arange(5, n_genes, 5),
            xmin=-0.5,
            xmax=n_strains - 0.5,
            colors="gray",
            linewidths=0.3,
            alpha=0.2,
            zorder=3,
        )
    if n_strains > 25:
        ax_heat.vlines(
            np.arange(25, n_strains, 25) - 0.5,
            ymin=0,
            ymax=n_genes,
            colors="gray",
            linewidths=0.3,
            alpha=0.2,
            zorder=3,
        )

    rect = Rectangle(
        (x0 - 0.5, y0),
        (x1 - x0) + 1,
        (y1 - y0) + 1,
        fill=False,
        edgecolor=st_to_color[HIGHLIGHT_ST],
        linewidth=1.5,
        zorder=10,
    )
    ax_heat.add_patch(rect)

    ax_heat.set_xticks([])
    ax_heat.set_yticks(np.arange(n_genes) + 0.5)
    ax_heat.set_yticklabels([gene_label(g) for g in gene_order_clust], fontsize=4)

    ax_st.imshow(
        strip_arr,
        aspect="auto",
        cmap=strip_cmap,
        interpolation="none",
        origin="upper",
        extent=[-0.5, n_strains - 0.5, 1, 0],
        vmin=0,
        vmax=n_other,
    )
    ax_st.set_xlim(ax_heat.get_xlim())
    ax_st.set_xticks([])
    ax_st.set_yticks([])
    for s in ax_st.spines.values():
        s.set_visible(False)

    if n_strains > 25:
        ax_tree.vlines(
            np.arange(25, n_strains, 25) - 0.5,
            ymin=0,
            ymax=max_depth,
            colors="gray",
            linewidths=0.3,
            alpha=0.2,
            zorder=0,
        )
    lc = LineCollection(segs, colors="black", linewidths=0.4, zorder=2)
    ax_tree.add_collection(lc)
    ax_tree.set_xlim(ax_heat.get_xlim())
    ax_tree.set_ylim(max_depth * 1.02, -max_depth * 0.02)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])
    for s in ax_tree.spines.values():
        s.set_visible(False)

    ax_legend.axis("off")
    handles = [
        Patch(facecolor=st_to_color[st], label=f"ST{st}") for st in highlight_sts
    ]
    ax_legend.legend(
        handles=handles,
        loc="center",
        fontsize=6,
        frameon=False,
        ncol=6,
        columnspacing=0.8,
        handlelength=1.2,
        handletextpad=0.4,
        title=f"ST (>={ST_MIN_STRAINS} strains)",
        title_fontsize=8,
    )

    freq_ord = full_freq.loc[gene_order_clust].values
    ax_freq.barh(
        np.arange(n_genes) + 0.5,
        freq_ord,
        height=1.0,
        color="lightgray",
        edgecolor="none",
    )
    ax_freq.set_ylim(ax_heat.get_ylim())
    ax_freq.set_yticks([])
    ax_freq.set_xlim(0, 1)
    ax_freq.set_xlabel("freq in full dataset", fontsize=9)
    ax_freq.spines[["top", "right"]].set_visible(False)

    save_fig(fig, args.out_heatmap)
    print(f"wrote {args.out_heatmap}")

    # Supplementary CSV
    def _gene_name(gid):
        name = id_to_name.get(gid)
        if pd.isna(name) or name in (None, ""):
            return ""
        return str(name)

    out_csv = pathlib.Path(args.out_pa_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["geneId", "gene_name"] + list(tip_order))
        w.writerow(["", ""] + [strain_st.get(s, "") for s in tip_order])
        sub = pa.loc[gene_order_clust, tip_order].astype(int)
        for gid, row in zip(gene_order_clust, sub.values):
            w.writerow([int(gid), _gene_name(gid)] + row.tolist())
    print(f"wrote {out_csv} ({sub.shape[0]} genes x {sub.shape[1]} strains)")


if __name__ == "__main__":
    main()
