#!/usr/bin/env python3

"""
Build GO count and presence/absence matrices from parsed eggNOG GO annotations.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build GO matrices from parsed eggNOG GO annotations."
    )
    parser.add_argument(
        "--go_long_tsv",
        required=True,
        help="Path to eggnog_go_long.tsv.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for GO matrix TSV files.",
    )
    return parser.parse_args()


def main() -> None:
    """Run GO matrix construction."""
    args = parse_args()
    go_long_tsv = Path(args.go_long_tsv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    go_df = pd.read_csv(
        filepath_or_buffer=go_long_tsv,
        sep="\t",
        dtype=str,
    )

    if go_df.empty:
        raise ValueError(f"GO long table is empty: {go_long_tsv}")

    go_df = go_df.drop_duplicates()
    go_df["value"] = 1

    count_matrix = (
        go_df.groupby(["isolate_id", "go_term"], as_index=False)["value"]
        .sum()
        .pivot(index="isolate_id", columns="go_term", values="value")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    count_matrix.to_csv(
        path_or_buf=out_dir / "go_count_matrix.tsv",
        sep="\t",
        index=False,
    )

    presence_df = count_matrix.copy()
    go_columns = [
        column_name
        for column_name in presence_df.columns
        if column_name != "isolate_id"
    ]
    presence_df[go_columns] = (presence_df[go_columns] > 0).astype(int)
    presence_df.to_csv(
        path_or_buf=out_dir / "go_presence_absence_matrix.tsv",
        sep="\t",
        index=False,
    )

    go_summary = (
        go_df.groupby("go_term", as_index=False)
        .agg(
            n_isolates=("isolate_id", "nunique"),
            n_proteins=("value", "sum"),
        )
        .sort_values(
            by=["n_isolates", "n_proteins", "go_term"],
            ascending=[False, False, True],
        )
    )
    go_summary.to_csv(
        path_or_buf=out_dir / "go_summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
