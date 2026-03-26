#!/usr/bin/env python3

"""
Parse antiSMASH per-isolate results into combined TSV outputs.

This script searches each antiSMASH isolate directory for region GenBank files,
extracts region-level cluster annotations, and writes long-format and matrix
outputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse antiSMASH result folders."
    )
    parser.add_argument(
        "--antismash_dir",
        required=True,
        help="Directory containing per-isolate antiSMASH result folders.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for parsed antiSMASH TSV files.",
    )
    return parser.parse_args()


def collect_region_files(isolate_dir: Path) -> List[Path]:
    """
    Collect antiSMASH region GenBank files for an isolate.

    Parameters
    ----------
    isolate_dir : Path
        Per-isolate antiSMASH result directory.

    Returns
    -------
    List[Path]
        Sorted region GenBank files.
    """
    patterns = [
        "*.region*.gbk",
        "*.region*.genbank",
        "*region*.gbk",
    ]

    region_files: List[Path] = []
    for pattern in patterns:
        region_files.extend(isolate_dir.glob(pattern))

    return sorted(set(region_files))


def get_feature_qualifier(feature, key: str) -> str:
    """
    Return a GenBank feature qualifier as a semicolon-joined string.

    Parameters
    ----------
    feature : SeqFeature
        Biopython feature object.
    key : str
        Qualifier key.

    Returns
    -------
    str
        Joined qualifier value.
    """
    values = feature.qualifiers.get(key, [])
    if not values:
        return ""
    return ";".join(str(value) for value in values)


def parse_region_file(isolate_id: str, region_file: Path) -> List[dict]:
    """
    Parse a single antiSMASH region GenBank file.

    Parameters
    ----------
    isolate_id : str
        Isolate identifier.
    region_file : Path
        Region GenBank file.

    Returns
    -------
    List[dict]
        Region-level annotation records.
    """
    records = []

    for record in SeqIO.parse(handle=str(region_file), format="genbank"):
        for feature in record.features:
            if feature.type != "region":
                continue

            product = get_feature_qualifier(feature=feature, key="product")
            category = get_feature_qualifier(feature=feature, key="category")
            contig_edge = get_feature_qualifier(feature=feature, key="contig_edge")
            rules = get_feature_qualifier(feature=feature, key="rules")
            region_number = get_feature_qualifier(
                feature=feature,
                key="region_number",
            )

            records.append(
                {
                    "isolate_id": isolate_id,
                    "region_file": str(region_file),
                    "record_id": record.id,
                    "record_description": record.description,
                    "region_number": region_number,
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                    "product": product,
                    "category": category,
                    "contig_edge": contig_edge,
                    "rules": rules,
                }
            )

    return records


def main() -> None:
    """Run the antiSMASH parsing workflow."""
    args = parse_args()
    antismash_dir = Path(args.antismash_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_records = []

    for isolate_dir in sorted([path for path in antismash_dir.iterdir() if path.is_dir()]):
        region_files = collect_region_files(isolate_dir=isolate_dir)
        for region_file in region_files:
            all_records.extend(
                parse_region_file(
                    isolate_id=isolate_dir.name,
                    region_file=region_file,
                )
            )

    if not all_records:
        raise FileNotFoundError(
            f"No antiSMASH region files found under {antismash_dir}"
        )

    region_df = pd.DataFrame.from_records(all_records).drop_duplicates()
    region_df.to_csv(
        path_or_buf=out_dir / "antismash_region_long.tsv",
        sep="\t",
        index=False,
    )

    summary_df = (
        region_df.groupby("isolate_id", as_index=False)
        .agg(
            n_regions=("region_file", "count"),
            n_unique_products=("product", "nunique"),
        )
    )
    summary_df.to_csv(
        path_or_buf=out_dir / "antismash_isolate_summary.tsv",
        sep="\t",
        index=False,
    )

    product_long = region_df.loc[
        region_df["product"].fillna("").astype(str).str.len() > 0,
        ["isolate_id", "region_number", "product"],
    ].drop_duplicates()

    count_matrix = (
        product_long.assign(value=1)
        .groupby(["isolate_id", "product"], as_index=False)["value"]
        .sum()
        .pivot(index="isolate_id", columns="product", values="value")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    count_matrix.to_csv(
        path_or_buf=out_dir / "antismash_product_count_matrix.tsv",
        sep="\t",
        index=False,
    )

    presence_df = count_matrix.copy()
    product_columns = [
        column_name
        for column_name in presence_df.columns
        if column_name != "isolate_id"
    ]
    presence_df[product_columns] = (presence_df[product_columns] > 0).astype(int)
    presence_df.to_csv(
        path_or_buf=out_dir / "antismash_product_presence_absence_matrix.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
