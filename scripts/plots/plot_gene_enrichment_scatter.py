"""Per-gene enrichment scatter: frequency in high-load vs other strains.

For each focal gene plot frequency in high-load strains (x) vs frequency in
the rest (y). Marker size encodes the effective number of STs the gene
participates in. A shaded wedge highlights `freq_high_load > ENRICHMENT_FOLD
* freq_other`. Curated focal genes are labeled with names from the highlight
CSV.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

from utils import save_fig


MIN_ST_SIZE = 10      # min isolates per ST for participation ratio
ENRICHMENT_FOLD = 2   # highlight wedge: freq_high_load > fold * freq_other


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pa", required=True, help="Focal-gene presence/absence CSV")
    p.add_argument("--gene-info", required=True, help="Focal gene_info CSV")
    p.add_argument("--metadata", required=True, help="Horesh F1 genome metadata CSV")
    p.add_argument("--highlight", required=True, help="Curated highlight CSV (id,name,color)")
    p.add_argument("--out-fig-base", required=True, help="Output figure path stem (no extension)")
    p.add_argument("--high-load-threshold", type=int, required=True, help="Gene-count threshold marking 'high-load' isolates")
    return p.parse_args()


def main():
    args = parse_args()

    pa = pd.read_csv(args.pa, index_col=0)
    pa.index = pa.index.astype(str)

    meta = pd.read_csv(args.metadata)
    meta = meta.dropna(subset=["name_in_presence_absence"])
    strain_meta = meta.set_index("name_in_presence_absence")

    gene_info = pd.read_csv(args.gene_info)
    gene_info["geneId"] = gene_info["geneId"].astype(str)
    gene_info = gene_info.set_index("geneId")

    highlight = pd.read_csv(args.highlight)
    highlight["id"] = highlight["id"].astype(str)

    # High-load / other strains
    total_count = pa.sum(axis=0)
    high_load_strains = total_count[total_count >= args.high_load_threshold].index
    other_strains = total_count[total_count < args.high_load_threshold].index
    print(f"High-load strains (>= {args.high_load_threshold} genes): {len(high_load_strains)}")
    print(f"Other strains: {len(other_strains)}")

    freq_high_load = pa[high_load_strains].mean(axis=1)
    freq_other = pa[other_strains].mean(axis=1)

    # ST participation ratio
    strains_in_pa = strain_meta.reindex(pa.columns).dropna(subset=["ST"]).copy()
    strains_in_pa["ST"] = strains_in_pa["ST"].astype(str).str.strip()
    strains_in_pa = strains_in_pa[~strains_in_pa["ST"].str.endswith("~")]
    st_counts = strains_in_pa["ST"].value_counts()
    valid_sts = st_counts[st_counts >= MIN_ST_SIZE].index
    print(f"STs with >= {MIN_ST_SIZE} isolates: {len(valid_sts)}")

    st_freq = pd.DataFrame(index=pa.index, columns=valid_sts, dtype=float)
    for st in valid_sts:
        st_strains = strains_in_pa[strains_in_pa["ST"] == st].index
        st_freq[st] = pa[st_strains].mean(axis=1)

    # Effective # of STs the gene "participates" in: (sum f)^2 / sum(f^2)
    sum_f = st_freq.sum(axis=1)
    sum_f2 = (st_freq**2).sum(axis=1)
    participation = (sum_f**2) / sum_f2.replace(0, np.nan)

    genes = pd.DataFrame(
        {
            "freq_high_load": freq_high_load,
            "freq_other": freq_other,
            "participation": participation,
        }
    )
    genes["gene_name"] = gene_info["gene_name"].reindex(genes.index)

    in_region = genes["freq_high_load"] > ENRICHMENT_FOLD * genes["freq_other"]
    print(f"Genes with freq_high_load > {ENRICHMENT_FOLD} * freq_other: {in_region.sum()}")
    print(genes.loc[in_region].sort_values("freq_high_load", ascending=False).head(30))

    # Plot
    size_scale = 5

    fig = plt.figure(figsize=(8, 7.2))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        wspace=0.05,
        hspace=0.05,
    )
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

    bins = np.linspace(0, 1, 51)
    ax_top.hist(
        genes["freq_high_load"].dropna(),
        bins=bins,
        color="lightgray",
        edgecolor="gray",
        linewidth=0.3,
    )
    ax_right.hist(
        genes["freq_other"].dropna(),
        bins=bins,
        orientation="horizontal",
        color="lightgray",
        edgecolor="gray",
        linewidth=0.3,
    )
    ax_top.set_ylabel("# genes")
    ax_right.set_xlabel("# genes")
    ax_top.tick_params(labelbottom=False)
    ax_right.tick_params(labelleft=False)
    ax_top.spines[["top", "right"]].set_visible(False)
    ax_right.spines[["top", "right"]].set_visible(False)

    # Enrichment wedge
    ax.fill(
        [0, 1, 1],
        [0, 0, 1.0 / ENRICHMENT_FOLD],
        color="orange",
        alpha=0.12,
        zorder=0,
    )
    ax.plot(
        [0, 1],
        [0, 1.0 / ENRICHMENT_FOLD],
        "k--",
        linewidth=0.8,
        alpha=0.6,
        zorder=1,
    )

    ax.scatter(
        genes["freq_high_load"],
        genes["freq_other"],
        s=genes["participation"] * size_scale,
        alpha=0.4,
        c="lightgray",
        edgecolor="gray",
        linewidth=0.3,
        zorder=2,
    )

    texts = []
    for _, row in highlight.iterrows():
        gid = row["id"]
        if gid not in genes.index:
            print(f"  warning: highlight id {gid} ({row['name']}) not in presence/absence")
            continue
        g = genes.loc[gid]
        ax.scatter(
            g["freq_high_load"],
            g["freq_other"],
            s=g["participation"] * size_scale,
            c=row["color"],
            edgecolor="black",
            linewidth=0.6,
            zorder=5,
        )
        texts.append(
            ax.text(
                g["freq_high_load"],
                g["freq_other"],
                row["name"],
                fontsize=7,
                zorder=6,
            )
        )

    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
        expand=(1.4, 1.6),
        force_text=(0.4, 0.6),
    )

    ax.set_xlabel("Frequency in high-load strains")
    ax.set_ylabel("Frequency in other strains")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)

    size_examples = [1, 5, 20, int(len(valid_sts))]
    size_handles = [
        plt.scatter(
            [],
            [],
            s=s * size_scale,
            color="lightgray",
            edgecolor="gray",
            linewidth=0.3,
            label=str(s),
        )
        for s in size_examples
    ]
    ax.legend(
        handles=size_handles,
        title="Effective # STs",
        loc="center left",
        labelspacing=1,
        borderpad=1.0,
        frameon=True,
        fontsize=7,
        title_fontsize=9,
    )

    plt.tight_layout()
    save_fig(fig, args.out_fig_base)
    print(f"Saved {args.out_fig_base}.(png|pdf)")


if __name__ == "__main__":
    main()
