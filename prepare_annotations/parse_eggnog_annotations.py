#!/usr/bin/env python3

"""
Parse eggNOG-mapper annotation files into combined TSV outputs.

This script reads all eggNOG-mapper annotation files from per-isolate result
directories, combines them into a single master table, and writes long-format
tables for KOs, GO terms, EC numbers, and COG categories.

Expected input layout
---------------------
<eggnog_results_dir>/
    <isolate_id>/
        <isolate_id>.emapper.annotations

Outputs
-------
- eggnog_master_annotations.tsv
- eggnog_ko_long.tsv
- eggnog_go_long.tsv
- eggnog_ec_long.tsv
- eggnog_cog_long.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse eggNOG-mapper annotation files."
    )
    parser.add_argument(
        "--eggnog_dir",
        required=True,
        help="Directory containing per-isolate eggNOG result folders.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for parsed eggNOG TSV files.",
    )
    return parser.parse_args()


def split_field(value: object, separator: str = ",") -> List[str]:
    """
    Split a delimited annotation field into a clean list of values.

    Parameters
    ----------
    value : object
        Raw field value.
    separator : str, optional
        Separator character, by default ",".

    Returns
    -------
    List[str]
        Clean list of non-empty values.
    """
    if pd.isna(value):
        return []

    text = str(value).strip()
    if text in {"", "-", "None", "NA", "nan"}:
        return []

    return [item.strip() for item in text.split(separator) if item.strip()]


def load_emapper_annotations(annotation_path: Path) -> pd.DataFrame:
    """
    Load a single eggNOG annotation file.

    Parameters
    ----------
    annotation_path : Path
        Path to a ``*.emapper.annotations`` file.

    Returns
    -------
    pd.DataFrame
        Parsed annotation table.
    """
    df = pd.read_csv(
        filepath_or_buffer=annotation_path,
        sep="\t",
        comment="#",
        dtype=str,
        low_memory=False,
    )
    return df


def collect_annotation_files(eggnog_dir: Path) -> List[Tuple[str, Path]]:
    """
    Collect per-isolate eggNOG annotation files.

    Parameters
    ----------
    eggnog_dir : Path
        Root eggNOG results directory.

    Returns
    -------
    List[Tuple[str, Path]]
        List of ``(isolate_id, annotation_file)`` tuples.
    """
    annotation_files: List[Tuple[str, Path]] = []

    for isolate_dir in sorted(
        [path for path in eggnog_dir.iterdir() if path.is_dir()]
    ):
        matches = list(isolate_dir.glob("*.emapper.annotations"))
        if not matches:
            continue
        annotation_files.append((isolate_dir.name, matches[0]))

    return annotation_files


def build_long_table(
    df: pd.DataFrame,
    source_column: str,
    value_name: str,
    separator: str = ",",
) -> pd.DataFrame:
    """
    Convert a delimited annotation column into a long-format table.

    Parameters
    ----------
    df : pd.DataFrame
        Master annotation table.
    source_column : str
        Column to split.
    value_name : str
        Name of the long-format value column.
    separator : str, optional
        Delimiter used in the source column, by default ",".

    Returns
    -------
    pd.DataFrame
        Long-format annotation table.
    """
    records = []

    if source_column not in df.columns:
        return pd.DataFrame(
            columns=[
                "isolate_id",
                "#query",
                value_name,
            ]
        )

    for _, row in df.iterrows():
        entries = split_field(value=row[source_column], separator=separator)
        for entry in entries:
            records.append(
                {
                    "isolate_id": row["isolate_id"],
                    "#query": row["#query"],
                    value_name: entry,
                }
            )

    if not records:
        return pd.DataFrame(
            columns=[
                "isolate_id",
                "#query",
                value_name,
            ]
        )

    return pd.DataFrame.from_records(records).drop_duplicates()


def build_cog_long_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expand the eggNOG COG category field into one row per category letter.

    Parameters
    ----------
    df : pd.DataFrame
        Master annotation table.

    Returns
    -------
    pd.DataFrame
        Long-format COG category table.
    """
    records = []

    if "COG_category" not in df.columns:
        return pd.DataFrame(
            columns=[
                "isolate_id",
                "#query",
                "cog_category",
            ]
        )

    for _, row in df.iterrows():
        value = row["COG_category"]
        if pd.isna(value):
            continue

        text = str(value).strip()
        if text in {"", "-", "NA", "nan"}:
            continue

        for category in text:
            if category.strip():
                records.append(
                    {
                        "isolate_id": row["isolate_id"],
                        "#query": row["#query"],
                        "cog_category": category,
                    }
                )

    if not records:
        return pd.DataFrame(
            columns=[
                "isolate_id",
                "#query",
                "cog_category",
            ]
        )

    return pd.DataFrame.from_records(records).drop_duplicates()


def main() -> None:
    """Run the eggNOG parsing workflow."""
    args = parse_args()
    eggnog_dir = Path(args.eggnog_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    annotation_files = collect_annotation_files(eggnog_dir=eggnog_dir)
    if not annotation_files:
        raise FileNotFoundError(
            f"No *.emapper.annotations files found under {eggnog_dir}"
        )

    combined_frames = []

    for isolate_id, annotation_path in annotation_files:
        df = load_emapper_annotations(annotation_path=annotation_path)
        df["isolate_id"] = isolate_id
        df["annotation_file"] = str(annotation_path)
        combined_frames.append(df)

    master_df = pd.concat(combined_frames, ignore_index=True)

    master_out = out_dir / "eggnog_master_annotations.tsv"
    master_df.to_csv(path_or_buf=master_out, sep="\t", index=False)

    ko_long = build_long_table(
        df=master_df,
        source_column="KEGG_ko",
        value_name="ko_id",
    )
    ko_long.to_csv(
        path_or_buf=out_dir / "eggnog_ko_long.tsv",
        sep="\t",
        index=False,
    )

    go_long = build_long_table(
        df=master_df,
        source_column="GOs",
        value_name="go_term",
    )
    go_long.to_csv(
        path_or_buf=out_dir / "eggnog_go_long.tsv",
        sep="\t",
        index=False,
    )

    ec_long = build_long_table(
        df=master_df,
        source_column="EC",
        value_name="ec_number",
    )
    ec_long.to_csv(
        path_or_buf=out_dir / "eggnog_ec_long.tsv",
        sep="\t",
        index=False,
    )

    cog_long = build_cog_long_table(df=master_df)
    cog_long.to_csv(
        path_or_buf=out_dir / "eggnog_cog_long.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
