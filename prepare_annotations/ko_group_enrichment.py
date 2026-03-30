#!/usr/bin/env python3

"""
Perform KO enrichment testing between two isolate groups.

This script compares KO presence/absence between two isolate groups using
Fisher's exact test and applies Benjamini-Hochberg multiple testing
correction.

Input
-----
- KO presence/absence matrix TSV
- Metadata TSV containing isolate IDs and a grouping column

Outputs
-------
- ko_group_enrichment.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Perform KO enrichment between two isolate groups."
    )
    parser.add_argument(
        "--ko_matrix_tsv",
        required=True,
        help="Path to KO presence/absence matrix TSV.",
    )
    parser.add_argument(
        "--metadata_tsv",
        required=True,
        help="Path to isolate metadata TSV.",
    )
    parser.add_argument(
        "--group_column",
        required=True,
        help="Metadata column containing group labels.",
    )
    parser.add_argument(
        "--group_a",
        required=True,
        help="First group label.",
    )
    parser.add_argument(
        "--group_b",
        required=True,
        help="Second group label.",
    )
    parser.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV path.",
    )
    return parser.parse_args()


def main() -> None:
    """
    Run KO enrichment testing.
    """
    args = parse_args()

    ko_df = pd.read_csv(args.ko_matrix_tsv, sep="\t", dtype=str)
    metadata_df = pd.read_csv(args.metadata_tsv, sep="\t", dtype=str)

    if "isolate_id" not in ko_df.columns:
        raise ValueError("KO matrix must contain an 'isolate_id' column.")

    if "isolate_id" not in metadata_df.columns:
        raise ValueError("Metadata TSV must contain an 'isolate_id' column.")

    if args.group_column not in metadata_df.columns:
        raise ValueError(
            f"Metadata TSV does not contain group column '{args.group_column}'."
        )

    ko_df = ko_df.set_index("isolate_id")
    ko_df = ko_df.apply(pd.to_numeric)

    metadata_df = metadata_df.loc[
        metadata_df[args.group_column].isin([args.group_a, args.group_b])
    ].copy()

    common_isolates = sorted(set(ko_df.index).intersection(metadata_df["isolate_id"]))
    if not common_isolates:
        raise ValueError("No shared isolate IDs between KO matrix and metadata.")

    ko_df = ko_df.loc[common_isolates]
    metadata_df = metadata_df.set_index("isolate_id").loc[common_isolates]

    group_a_isolates = metadata_df.index[
        metadata_df[args.group_column] == args.group_a
    ].tolist()
    group_b_isolates = metadata_df.index[
        metadata_df[args.group_column] == args.group_b
    ].tolist()

    if not group_a_isolates or not group_b_isolates:
        raise ValueError("One or both groups contain no isolates.")

    results = []

    for ko_id in ko_df.columns:
        a_present = int((ko_df.loc[group_a_isolates, ko_id] > 0).sum())
        a_absent = len(group_a_isolates) - a_present
        b_present = int((ko_df.loc[group_b_isolates, ko_id] > 0).sum())
        b_absent = len(group_b_isolates) - b_present

        table = np.array(
            [
                [a_present, a_absent],
                [b_present, b_absent],
            ]
        )

        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")

        prevalence_a = a_present / len(group_a_isolates)
        prevalence_b = b_present / len(group_b_isolates)

        results.append(
            {
                "ko_id": ko_id,
                "group_a": args.group_a,
                "group_b": args.group_b,
                "group_a_n": len(group_a_isolates),
                "group_b_n": len(group_b_isolates),
                "group_a_present": a_present,
                "group_a_absent": a_absent,
                "group_b_present": b_present,
                "group_b_absent": b_absent,
                "group_a_prevalence": prevalence_a,
                "group_b_prevalence": prevalence_b,
                "prevalence_difference": prevalence_a - prevalence_b,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
            }
        )

    results_df = pd.DataFrame(results)

    reject, p_adj, _, _ = multipletests(
        pvals=results_df["p_value"].fillna(1.0).values,
        method="fdr_bh",
    )
    results_df["p_adj_bh"] = p_adj
    results_df["significant_bh"] = reject

    results_df = results_df.sort_values(
        by=["p_adj_bh", "p_value", "prevalence_difference"],
        ascending=[True, True, False],
    )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
