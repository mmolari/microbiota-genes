#!/usr/bin/env python3
"""
Filter homology hits (blastn) and assign each focal-gene consensus to
its best Horesh reference(s).

Input contract
--------------
A 12-column tabular file with the same columns as `blastn -outfmt 6`:
    qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore

BLAST emits one row per HSP, so a single (query, subject) hit can span
multiple rows. We aggregate per (qseqid, sseqid) before applying coverage
thresholds: qcov / scov are computed as the *union* of HSP intervals on the
query and subject respectively, divided by qlen / slen. This matches BLAST's
own `qcovs` definition on the query side and applies the same logic to the
subject side.

Coverage note
-------------
Focal-gene consensus sequences (built from MSAs) are often longer than the
individual pangenome reference genes, and can split into multiple HSPs against
a single subject. Aggregating intervals before thresholding avoids spuriously
dropping multi-HSP hits where each individual HSP is below threshold but the
union of HSPs is not.

Output
------
mapping.csv   - one row per (query, ref) pair passing all thresholds.
                Columns:
                  query_id, ref_id, pident, qcov, scov,
                  bitscore, evalue, n_hsps, hit_rank
                pident: length-weighted mean across HSPs;
                bitscore: sum across HSPs; evalue: min across HSPs;
                hit_rank: 1 = best ref for this query (by aggregated bitscore).
unmatched.txt - one query ID per line for queries with no passing hit.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
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

OUT_COLS = [
    "query_id",
    "ref_id",
    "pident",
    "qcov",
    "scov",
    "bitscore",
    "evalue",
    "n_hsps",
    "hit_rank",
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


def union_length(intervals):
    """Length of the union of 1-based inclusive intervals; orientation-agnostic."""
    norm = sorted((min(s, e), max(s, e)) for s, e in intervals)
    total, cur_end = 0, 0
    for s, e in norm:
        s = max(s, cur_end + 1)
        if s <= e:
            total += e - s + 1
            cur_end = e
    return total


def aggregate_hits(df):
    """Collapse per-HSP rows into one row per (qseqid, sseqid)."""
    records = []
    for (qid, sid), grp in df.groupby(["qseqid", "sseqid"], sort=False):
        qlen = grp["qlen"].iloc[0]
        slen = grp["slen"].iloc[0]
        q_union = union_length(zip(grp["qstart"], grp["qend"]))
        s_union = union_length(zip(grp["sstart"], grp["send"]))
        records.append(
            {
                "qseqid": qid,
                "sseqid": sid,
                "pident": float(np.average(grp["pident"], weights=grp["length"])),
                "qcov": q_union / qlen * 100,
                "scov": s_union / slen * 100,
                "bitscore": float(grp["bitscore"].sum()),
                "evalue": float(grp["evalue"].min()),
                "n_hsps": int(len(grp)),
            }
        )
    return pd.DataFrame.from_records(records)


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
        agg = pd.DataFrame(
            columns=["qseqid", "sseqid", "pident", "qcov", "scov",
                     "bitscore", "evalue", "n_hsps"]
        )
        passing = agg.copy()
    else:
        # --- aggregate per (qseqid, sseqid) before thresholding ---
        agg = aggregate_hits(df)
        agg["pident"] = agg["pident"].round(2)
        agg["qcov"] = agg["qcov"].round(2)
        agg["scov"] = agg["scov"].round(2)

        passing = agg[
            (agg["pident"] >= args.min_pident)
            & (agg["qcov"] >= args.min_qcovs)
            & (agg["scov"] >= args.min_scovs)
        ].copy()

    print(
        f"HSPs: {len(df)}; aggregated hits: {len(agg)}; passing "
        f"pident>={args.min_pident} qcov>={args.min_qcovs} scov>={args.min_scovs}: "
        f"{len(passing)}",
        file=sys.stderr,
    )

    # --- rank hits within each query (best aggregated bitscore = rank 1) ---
    if not passing.empty:
        passing = passing.sort_values(["qseqid", "bitscore"], ascending=[True, False])
        passing["hit_rank"] = passing.groupby("qseqid").cumcount() + 1
    else:
        passing["hit_rank"] = pd.Series(dtype=int)

    # --- build final output ---
    out = passing.rename(columns={"qseqid": "query_id", "sseqid": "ref_id"})
    out = out.reindex(columns=OUT_COLS)

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

    # --- summaries of multi-hit / multi-HSP cases ---
    if not passing.empty:
        multi = passing[passing["hit_rank"] > 1]["qseqid"].nunique()
        if multi:
            print(
                f"  {multi} queries have >1 passing ref (possible paralogs or multi-gene consensus)",
                file=sys.stderr,
            )
        multi_hsp = int((passing["n_hsps"] > 1).sum())
        if multi_hsp:
            print(
                f"  {multi_hsp} passing (query, ref) pairs aggregated >1 HSP",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
