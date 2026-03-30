#!/usr/bin/env python3

"""
Plot a PGPR presence/absence heatmap.

This script reads a PGPR presence/absence matrix TSV and generates a heatmap
PNG plus an ordered TSV for downstream review.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Plot a PGPR presence/absence heatmap."
    )
    parser.add_argument(
        "--pgpr_matrix_tsv",
        required=True,
        help="Path to PGPR presence/absence matrix TSV.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for heatmap outputs.",
    )
    return parser.parse_args()


def main() -> None:
    """
    Run the PGPR heatmap workflow.
    """
    args = parse_args()

    pgpr_df = pd.read_csv(args.pgpr_matrix_tsv, sep="\t", dtype=str)
    if "isolate_id" not in pgpr_df.columns:
        raise ValueError("PGPR matrix must contain an 'isolate_id' column.")

    pgpr_df = pgpr_df.set_index("isolate_id")
    pgpr_df = pgpr_df.apply(pd.to_numeric)

    row_sums = pgpr_df.sum(axis=1)
    pgpr_df = pgpr_df.loc[row_sums.sort_values(ascending=False).index]

    col_sums = pgpr_df.sum(axis=0)
    pgpr_df = pgpr_df.loc[:, col_sums.sort_values(ascending=False).index]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ordered_tsv = out_dir / "pgpr_presence_absence_matrix_ordered.tsv"
    pgpr_df.reset_index().to_csv(ordered_tsv, sep="\t", index=False)

    fig_width = max(6, 1 + 0.8 * pgpr_df.shape[1])
    fig_height = max(4, 1 + 0.5 * pgpr_df.shape[0])

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    im = ax.imshow(pgpr_df.values, aspect="auto")

    ax.set_xticks(range(pgpr_df.shape[1]))
    ax.set_xticklabels(pgpr_df.columns, rotation=45, ha="right")
    ax.set_yticks(range(pgpr_df.shape[0]))
    ax.set_yticklabels(pgpr_df.index)

    ax.set_xlabel("PGPR gene group")
    ax.set_ylabel("Isolate")
    ax.set_title("PGPR gene-group presence/absence heatmap")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Presence / absence")

    fig.tight_layout()
    fig.savefig(out_dir / "pgpr_presence_absence_heatmap.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
