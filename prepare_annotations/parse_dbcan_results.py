#!/usr/bin/env python3

"""
Parse dbCAN result folders into combined TSV outputs.

This script reads dbCAN ``overview.txt`` files from per-isolate result
directories, combines them into a single master table, extracts CAZy family
assignments, and writes count and presence/absence matrices.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Set

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse dbCAN CAZyme annotation results."
    )
    parser.add_argument(
        "--dbcan_dir",
        required=True,
        help="Directory containing per-isolate dbCAN result folders.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for parsed dbCAN TSV files.",
    )
    return parser.parse_args()


def detect_overview_file(isolate_dir: Path) -> Path | None:
    """
    Detect the dbCAN overview file for an isolate.

    Parameters
    ----------
    isolate_dir : Path
        Per-isolate dbCAN result directory.

    Returns
    -------
    Path | None
        Path to the overview file, or None if not found.
    """
    candidates = [
        isolate_dir / "overview.txt",
        isolate_dir / "overview.tsv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def load_dbcan_overview(overview_path: Path) -> pd.DataFrame:
    """
    Load a dbCAN overview file.

    Parameters
    ----------
    overview_path : Path
        Path to dbCAN overview file.

    Returns
    -------
    pd.DataFrame
        Parsed overview table.
    """
    return pd.read_csv(
        filepath_or_buffer=overview_path,
        sep="\t",
        dtype=str,
        low_memory=False,
    )


def extract_cazy_families(text: object) -> List[str]:
    """
    Extract CAZy family identifiers from a dbCAN annotation field.

    Parameters
    ----------
    text : object
        Raw annotation field.

    Returns
    -------
    List[str]
        Sorted unique CAZy family identifiers.
    """
    if pd.isna(text):
        return []

    value = str(text)
    matches = re.findall(
        pattern=r"\b(?:GH|GT|CBM|CE|PL|AA)\d+(?:_\d+)?\b",
        string=value,
    )
    return sorted(set(matches))


def collect_dbcan_tables(dbcan_dir: Path) -> pd.DataFrame:
    """
    Collect and combine dbCAN overview tables.

    Parameters
    ----------
    dbcan_dir : Path
        Root dbCAN results directory.

    Returns
    -------
    pd.DataFrame
        Combined dbCAN overview table.
    """
    frames = []

    for isolate_dir in sorted([path for path in dbcan_dir.iterdir() if path.is_dir()]):
        overview_path = detect_overview_file(isolate_dir=isolate_dir)
        if overview_path is None:
            continue

        df = load_dbcan_overview(overview_path=overview_path)
        df["isolate_id"] = isolate_dir.name
        df["overview_file"] = str(overview_path)
        frames.append(df)

    if not frames:
        raise FileNotFoundError(
            f"No dbCAN overview files found under {dbcan_dir}"
        )

    return pd.concat(frames, ignore_index=True)


def pick_gene_column(df: pd.DataFrame) -> str:
    """
    Pick the most likely gene identifier column from a dbCAN table.

    Parameters
    ----------
    df : pd.DataFrame
        Combined dbCAN table.

    Returns
    -------
    str
        Chosen gene identifier column name.
    """
    candidates = ["Gene ID", "GeneID", "Protein ID", "Protein_ID"]
    for candidate in candidates:
        if candidate in df.columns:
            return candidate

    return df.columns[0]


def main() -> None:
    """Run the dbCAN parsing workflow."""
    args = parse_args()
    dbcan_dir = Path(args.dbcan_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    master_df = collect_dbcan_tables(dbcan_dir=dbcan_dir)
    gene_column = pick_gene_column(df=master_df)

    master_df.to_csv(
        path_or_buf=out_dir / "dbcan_master_annotations.tsv",
        sep="\t",
        index=False,
    )

    family_records = []

    source_columns = [
        column_name
        for column_name in ["HMMER", "dbCAN_sub", "DIAMOND"]
        if column_name in master_df.columns
    ]

    for _, row in master_df.iterrows():
        families: Set[str] = set()
        for source_column in source_columns:
            families.update(extract_cazy_families(text=row[source_column]))

        for family in sorted(families):
            family_records.append(
                {
                    "isolate_id": row["isolate_id"],
                    "gene_id": row[gene_column],
                    "cazy_family": family,
                }
            )

    family_df = pd.DataFrame.from_records(family_records).drop_duplicates()
    family_df.to_csv(
        path_or_buf=out_dir / "dbcan_cazy_family_long.tsv",
        sep="\t",
        index=False,
    )

    count_matrix = (
        family_df.assign(value=1)
        .groupby(["isolate_id", "cazy_family"], as_index=False)["value"]
        .sum()
        .pivot(index="isolate_id", columns="cazy_family", values="value")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    count_matrix.to_csv(
        path_or_buf=out_dir / "dbcan_cazy_family_count_matrix.tsv",
        sep="\t",
        index=False,
    )

    presence_df = count_matrix.copy()
    family_columns = [
        column_name
        for column_name in presence_df.columns
        if column_name != "isolate_id"
    ]
    presence_df[family_columns] = (presence_df[family_columns] > 0).astype(int)
    presence_df.to_csv(
        path_or_buf=out_dir / "dbcan_cazy_family_presence_absence_matrix.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
