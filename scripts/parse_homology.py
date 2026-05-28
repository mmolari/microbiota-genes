#!/usr/bin/env python3
"""
Filter homology hits (blastn) and assign each focal-gene consensus to
its best Horesh reference(s).

Input contract
--------------
A 12-column tabular file with the same columns as `blastn -outfmt 6`:
    qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
For mmseqs2, use `easy-search --format-output "query,target,pident,alnlen,qlen,
tlen,qstart,qend,tstart,tend,evalue,bits"`. With `--search-type 2`, `pident`,
`length`, `qlen` and `slen` are all in amino-acid units, so the qcov/scov
arithmetic below remains unit-consistent.

Coverage note
-------------
Focal-gene consensus sequences (built from MSAs) are often longer than the individual
pangenome reference genes.  Because of this, *subject coverage* (scovs = how much of
the reference gene is covered) is the more meaningful filter for a "complete match",
whereas query coverage (qcovs) flags whether the whole consensus maps to one reference.
Both thresholds are applied independently.

Output
------
mapping.csv   - all hits passing all thresholds, one row per (query, ref) pair.
                Columns:
                  query_id, ref_id, pident, qcov, scov,
                  qstart, qend, sstart, send, evalue, bitscore,
                  hit_rank          # 1 = best hit for this query
unmatched.txt - one query ID per line for queries with no passing hit
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


BLAST_COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hits", help="Homology hits TSV in blast -outfmt 6 column order")
    p.add_argument("--query-fa", help="Query FASTA (used to get full query ID list)")
    p.add_argument("--mapping", help="Output CSV with passing hits")
    p.add_argument("--unmatched", help="Output text file with unmatched query IDs")
    p.add_argument("--min-pident", type=float, help="Min percent identity")
    p.add_argument("--min-qcovs", type=float, help="Min query coverage")
    p.add_argument("--min-scovs", type=float, help="Min subject coverage")
    return p.parse_args()


def read_query_ids(fa_path):
    ids = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    args = parse_args()

    # --- load all query IDs so we can report unmatched ones ---
    all_query_ids = read_query_ids(args.query_fa)

    # --- load hits ---
    hits_path = Path(args.hits)
    if hits_path.stat().st_size == 0:
        df = pd.DataFrame(columns=BLAST_COLS)
    else:
        df = pd.read_csv(
            hits_path,
            sep="\t",
            header=None,
            names=BLAST_COLS,
            dtype={"qseqid": str, "sseqid": str},
        )

    if df.empty:
        print(
            "WARNING: hits file is empty — all queries will be unmatched.",
            file=sys.stderr,
        )
        passing = df.copy()
    else:
        # --- compute per-alignment coverage ---
        df["qcov"] = (df["length"] / df["qlen"] * 100).round(2)
        df["scov"] = (df["length"] / df["slen"] * 100).round(2)

        # --- apply thresholds ---
        passing = df[
            (df["pident"] >= args.min_pident)
            & (df["qcov"] >= args.min_qcovs)
            & (df["scov"] >= args.min_scovs)
        ].copy()

    print(
        f"Hits: {len(df)} total, {len(passing)} pass "
        f"pident>={args.min_pident}  qcov>={args.min_qcovs}  scov>={args.min_scovs}",
        file=sys.stderr,
    )

    # --- rank hits within each query (best bitscore = rank 1) ---
    if not passing.empty:
        passing = passing.sort_values(["qseqid", "bitscore"], ascending=[True, False])
        passing["hit_rank"] = passing.groupby("qseqid").cumcount() + 1

    # --- build final output ---
    out_cols = [
        "qseqid",
        "sseqid",
        "pident",
        "qcov",
        "scov",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "hit_rank",
    ]
    out = passing[out_cols].rename(columns={"qseqid": "query_id", "sseqid": "ref_id"})

    Path(args.mapping).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.mapping, index=False)

    # --- unmatched queries ---
    matched = set(passing["qseqid"].unique()) if not passing.empty else set()
    unmatched = [qid for qid in all_query_ids if qid not in matched]
    with open(args.unmatched, "w") as fh:
        for qid in unmatched:
            print(qid, file=fh)

    n_matched = len(all_query_ids) - len(unmatched)
    print(
        f"Queries: {len(all_query_ids)} total, "
        f"{n_matched} matched ({len(unmatched)} unmatched)",
        file=sys.stderr,
    )

    # --- summary of multi-hit queries ---
    if not passing.empty:
        multi = passing[passing["hit_rank"] > 1]["qseqid"].nunique()
        if multi:
            print(
                f"  {multi} queries have >1 passing hit (possible paralogs or multi-gene consensus)",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
