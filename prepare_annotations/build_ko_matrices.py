#!/usr/bin/env python3

"""
Build KO count and presence/absence matrices from parsed eggNOG KO annotations.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build KO matrices from parsed eggNOG KO annotations."
    )
    parser.add_argument(
        "--ko_long_tsv",
        required=True,
        help="Path to eggnog_ko_long.tsv.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for KO matrix TSV files.",
    )
    return parser.parse_args()


def main() -> None:
    """Run KO matrix construction."""
    args = parse_args()
    ko_long_tsv = Path(args.ko_long_tsv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ko_df = pd.read_csv(
        filepath_or_buffer=ko_long_tsv,
        sep="\t",
        dtype=str,
    )

    if ko_df.empty:
        raise ValueError(f"KO long table is empty: {ko_long_tsv}")

    ko_df = ko_df.drop_duplicates()
    ko_df["value"] = 1

    count_matrix = (
        ko_df.groupby(["isolate_id", "ko_id"], as_index=False)["value"]
        .sum()
        .pivot(index="isolate_id", columns="ko_id", values="value")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    count_matrix.to_csv(
        path_or_buf=out_dir / "ko_count_matrix.tsv",
        sep="\t",
        index=False,
    )

    presence_df = count_matrix.copy()
    ko_columns = [
        column_name
        for column_name in presence_df.columns
        if column_name != "isolate_id"
    ]
    presence_df[ko_columns] = (presence_df[ko_columns] > 0).astype(int)
    presence_df.to_csv(
        path_or_buf=out_dir / "ko_presence_absence_matrix.tsv",
        sep="\t",
        index=False,
    )

    ko_summary = (
        ko_df.groupby("ko_id", as_index=False)
        .agg(
            n_isolates=("isolate_id", "nunique"),
            n_proteins=("value", "sum"),
        )
        .sort_values(
            by=["n_isolates", "n_proteins", "ko_id"],
            ascending=[False, False, True],
        )
    )
    ko_summary.to_csv(
        path_or_buf=out_dir / "ko_summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
