"""Extract a focal-gene presence-absence sub-matrix from the Horesh PA matrix.

For each focal-gene query, presence is defined as the OR (max) across all
pangenome clusters assigned to that query via BLAST mapping.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract focal-gene presence-absence sub-matrix from Horesh PA matrix."
    )
    parser.add_argument("--mapping", help="Path to focal_vs_horesh_mapping.csv")
    parser.add_argument("--presence-absence", help="Path to full PA matrix CSV")
    parser.add_argument("--output", help="Path to output CSV")
    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    print("Loading BLAST mapping...", file=sys.stderr)
    mapping = pd.read_csv(args.mapping, usecols=["query_id", "ref_id"], dtype=str)
    ref_ids = set(mapping["ref_id"])
    print(
        f"  {len(mapping)} mapping rows, {len(ref_ids)} unique ref_ids, "
        f"{mapping['query_id'].nunique()} unique query_ids",
        file=sys.stderr,
    )

    print(
        "Loading presence-absence matrix (filtering to matched rows)...",
        file=sys.stderr,
    )
    pa = pd.read_csv(args.presence_absence, index_col=0)
    print(f"  Full PA matrix shape: {pa.shape}", file=sys.stderr)

    pa_sub = pa.loc[pa.index.isin(ref_ids)].copy()
    print(f"  Filtered PA matrix shape: {pa_sub.shape}", file=sys.stderr)
    del pa  # free memory

    print("Aggregating by query_id (OR across matched clusters)...", file=sys.stderr)
    rows = {}
    for query_id, group in mapping.groupby("query_id"):
        clusters = group["ref_id"].tolist()
        present_in_pa = [c for c in clusters if c in pa_sub.index]
        if not present_in_pa:
            continue
        rows[query_id] = pa_sub.loc[present_in_pa].max(axis=0)

    result = pd.DataFrame(rows).T
    result.index.name = "query_id"
    print(f"Result shape: {result.shape}", file=sys.stderr)

    print(f"Writing output to {args.output}...", file=sys.stderr)
    result.to_csv(args.output)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
