#!/usr/bin/env python3

"""
Screen eggNOG annotations for plant growth promotion-related gene groups.

This script searches eggNOG annotations for the requested groups:
pho, ding, pqq, nif, and hcn.

It uses preferred gene names, annotation descriptions, and KEGG KO strings.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


PGPR_PATTERNS: Dict[str, List[str]] = {
    "pho": [
        r"\bpho[a-z0-9_]*\b",
        r"\bpst[a-z0-9_]*\b",
        r"\bphn[a-z0-9_]*\b",
        r"phosphate",
        r"alkaline phosphatase",
        r"phosphate regulon",
    ],
    "ding": [
        r"\bding\b",
        r"ding protein",
    ],
    "pqq": [
        r"\bpqq[a-z0-9_]*\b",
        r"pyrroloquinoline quinone",
        r"glucose dehydrogenase",
    ],
    "nif": [
        r"\bnif[a-z0-9_]*\b",
        r"\bfix[a-z0-9_]*\b",
        r"nitrogen fixation",
        r"nitrogenase",
    ],
    "hcn": [
        r"\bhcn[a-z0-9_]*\b",
        r"hydrogen cyanide",
        r"cyanide synthase",
    ],
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Screen eggNOG annotations for PGPR gene groups."
    )
    parser.add_argument(
        "--eggnog_master_tsv",
        required=True,
        help="Path to eggnog_master_annotations.tsv.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for PGPR screen TSV files.",
    )
    return parser.parse_args()


def build_search_text(row: pd.Series) -> str:
    """
    Build a searchable text string from a master eggNOG annotation row.

    Parameters
    ----------
    row : pd.Series
        One annotation row.

    Returns
    -------
    str
        Combined search text.
    """
    candidate_columns = [
        "#query",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "PFAMs",
    ]

    values = []
    for column_name in candidate_columns:
        if column_name in row.index and pd.notna(row[column_name]):
            values.append(str(row[column_name]))

    return " | ".join(values)


def find_pgpr_hits(master_df: pd.DataFrame) -> pd.DataFrame:
    """
    Find PGPR-related hits in the eggNOG master annotation table.

    Parameters
    ----------
    master_df : pd.DataFrame
        Master eggNOG annotation table.

    Returns
    -------
    pd.DataFrame
        Long-format PGPR hit table.
    """
    records = []

    for _, row in master_df.iterrows():
        search_text = build_search_text(row=row).lower()

        for group_name, patterns in PGPR_PATTERNS.items():
            matched_patterns = [
                pattern
                for pattern in patterns
                if re.search(pattern=pattern, string=search_text, flags=re.IGNORECASE)
            ]

            if matched_patterns:
                records.append(
                    {
                        "isolate_id": row["isolate_id"],
                        "#query": row["#query"],
                        "pgpr_group": group_name,
                        "matched_patterns": ";".join(matched_patterns),
                        "preferred_name": row.get("Preferred_name", ""),
                        "description": row.get("Description", ""),
                        "kegg_ko": row.get("KEGG_ko", ""),
                    }
                )

    if not records:
        return pd.DataFrame(
            columns=[
                "isolate_id",
                "#query",
                "pgpr_group",
                "matched_patterns",
                "preferred_name",
                "description",
                "kegg_ko",
            ]
        )

    return pd.DataFrame.from_records(records).drop_duplicates()


def main() -> None:
    """Run the PGPR screen workflow."""
    args = parse_args()
    eggnog_master_tsv = Path(args.eggnog_master_tsv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    master_df = pd.read_csv(
        filepath_or_buffer=eggnog_master_tsv,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    hits_df = find_pgpr_hits(master_df=master_df)
    hits_df.to_csv(
        path_or_buf=out_dir / "pgpr_hits_long.tsv",
        sep="\t",
        index=False,
    )

    if hits_df.empty:
        empty_matrix = pd.DataFrame(columns=["isolate_id"])
        empty_matrix.to_csv(
            path_or_buf=out_dir / "pgpr_presence_absence_matrix.tsv",
            sep="\t",
            index=False,
        )
        return

    matrix_df = (
        hits_df.assign(value=1)
        .groupby(["isolate_id", "pgpr_group"], as_index=False)["value"]
        .max()
        .pivot(index="isolate_id", columns="pgpr_group", values="value")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    matrix_df.to_csv(
        path_or_buf=out_dir / "pgpr_presence_absence_matrix.tsv",
        sep="\t",
        index=False,
    )

    summary_df = (
        hits_df.groupby("isolate_id", as_index=False)
        .agg(
            n_pgpr_groups=("pgpr_group", "nunique"),
            n_pgpr_hits=("#query", "count"),
        )
        .sort_values(by=["n_pgpr_groups", "n_pgpr_hits"], ascending=[False, False])
    )
    summary_df.to_csv(
        path_or_buf=out_dir / "pgpr_isolate_summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
