#!/usr/bin/env python3

"""
Generate KO-based PCoA outputs and KO intersection summaries.

This script reads a KO presence/absence matrix, computes a Jaccard distance
matrix between isolates, performs classical multidimensional scaling (PCoA-like
ordination), and writes ordination tables and plots.

It also summarises KO intersections across isolates by identifying the isolate
combination in which each KO is present, then writing a long table and a
summary table of the most frequent isolate combinations.

Inputs
------
- KO presence/absence matrix TSV produced by build_ko_matrices.py

Outputs
-------
- ko_jaccard_distance_matrix.tsv
- ko_pcoa_coordinates.tsv
- ko_pcoa_variance_explained.tsv
- ko_pcoa_plot.png
- ko_isolate_membership_long.tsv
- ko_intersection_summary.tsv
- ko_top_intersections.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate KO PCoA outputs and KO intersection summaries."
    )
    parser.add_argument(
        "--ko_matrix_tsv",
        required=True,
        help="Path to KO presence/absence matrix TSV.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for PCoA and intersection outputs.",
    )
    parser.add_argument(
        "--top_n_intersections",
        type=int,
        default=20,
        help="Number of top KO intersections to plot.",
    )
    return parser.parse_args()


def load_ko_matrix(ko_matrix_tsv: Path) -> pd.DataFrame:
    """
    Load the KO presence/absence matrix.

    Parameters
    ----------
    ko_matrix_tsv : Path
        Path to the KO presence/absence matrix TSV.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by isolate_id with KO columns.
    """
    df = pd.read_csv(
        filepath_or_buffer=ko_matrix_tsv,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    if "isolate_id" not in df.columns:
        raise ValueError(
            "Expected 'isolate_id' column in KO matrix."
        )

    ko_df = df.copy()
    ko_df = ko_df.set_index("isolate_id")
    ko_df = ko_df.apply(pd.to_numeric)
    return ko_df


def compute_jaccard_distance_matrix(
    ko_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute a Jaccard distance matrix between isolates.

    Parameters
    ----------
    ko_df : pd.DataFrame
        Binary KO presence/absence matrix indexed by isolate.

    Returns
    -------
    pd.DataFrame
        Symmetric Jaccard distance matrix.
    """
    distances = pdist(X=ko_df.values, metric="jaccard")
    distance_matrix = squareform(X=distances)

    return pd.DataFrame(
        data=distance_matrix,
        index=ko_df.index,
        columns=ko_df.index,
    )


def classical_mds(
    distance_matrix: np.ndarray,
    n_components: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform classical multidimensional scaling.

    Parameters
    ----------
    distance_matrix : np.ndarray
        Square distance matrix.
    n_components : int, optional
        Number of dimensions to return, by default 2.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Coordinate matrix and eigenvalues.
    """
    n_samples = distance_matrix.shape[0]
    squared = distance_matrix ** 2

    centring = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    gram = -0.5 * centring @ squared @ centring

    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive = eigenvalues > 0
    eigenvalues = eigenvalues[positive]
    eigenvectors = eigenvectors[:, positive]

    if eigenvalues.size == 0:
        raise ValueError("No positive eigenvalues found in MDS decomposition.")

    coords = eigenvectors[:, :n_components] * np.sqrt(eigenvalues[:n_components])
    return coords, eigenvalues


def build_intersection_tables(ko_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build KO isolate-membership and intersection summary tables.

    Parameters
    ----------
    ko_df : pd.DataFrame
        Binary KO presence/absence matrix indexed by isolate.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Long membership table and intersection summary table.
    """
    membership_records = []
    summary_records = []

    for ko_id in ko_df.columns:
        present_isolates = ko_df.index[ko_df[ko_id] > 0].tolist()
        membership_key = "|".join(sorted(present_isolates))
        n_isolates = len(present_isolates)

        membership_records.append(
            {
                "ko_id": ko_id,
                "n_isolates_present": n_isolates,
                "isolate_membership_key": membership_key,
                "isolates_present": ";".join(sorted(present_isolates)),
            }
        )

    membership_df = pd.DataFrame.from_records(membership_records)

    summary_df = (
        membership_df.groupby(
            ["isolate_membership_key", "isolates_present", "n_isolates_present"],
            as_index=False,
        )
        .agg(n_kos=("ko_id", "count"))
        .sort_values(
            by=["n_kos", "n_isolates_present", "isolates_present"],
            ascending=[False, False, True],
        )
    )

    return membership_df, summary_df


def save_pcoa_plot(
    coords_df: pd.DataFrame,
    variance_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Save a PCoA scatter plot.

    Parameters
    ----------
    coords_df : pd.DataFrame
        PCoA coordinates table.
    variance_df : pd.DataFrame
        Variance explained table.
    output_path : Path
        Output PNG path.
    """
    pc1_var = variance_df.loc[0, "variance_explained_pct"]
    pc2_var = variance_df.loc[1, "variance_explained_pct"]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(coords_df["PC1"], coords_df["PC2"])

    for _, row in coords_df.iterrows():
        ax.text(
            x=row["PC1"],
            y=row["PC2"],
            s=row["isolate_id"],
            fontsize=8,
        )

    ax.set_xlabel(f"PCoA1 ({pc1_var:.2f}%)")
    ax.set_ylabel(f"PCoA2 ({pc2_var:.2f}%)")
    ax.set_title("PCoA of isolates by KO presence/absence")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def save_intersection_plot(
    summary_df: pd.DataFrame,
    output_path: Path,
    top_n: int,
) -> None:
    """
    Save a bar plot of the top KO intersection combinations.

    Parameters
    ----------
    summary_df : pd.DataFrame
        KO intersection summary table.
    output_path : Path
        Output PNG path.
    top_n : int
        Number of top intersections to plot.
    """
    top_df = summary_df.head(top_n).copy()

    if top_df.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No intersection data available", ha="center", va="center")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
        return

    top_df["plot_label"] = top_df["isolates_present"].str.replace(";", " | ")

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(top_df)), top_df["n_kos"])
    ax.set_xticks(range(len(top_df)))
    ax.set_xticklabels(top_df["plot_label"], rotation=90)
    ax.set_ylabel("Number of KOs")
    ax.set_title("Top KO isolate intersections")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    """
    Run the KO PCoA and intersection workflow.
    """
    args = parse_args()
    ko_matrix_tsv = Path(args.ko_matrix_tsv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ko_df = load_ko_matrix(ko_matrix_tsv=ko_matrix_tsv)

    distance_df = compute_jaccard_distance_matrix(ko_df=ko_df)
    distance_df.to_csv(
        path_or_buf=out_dir / "ko_jaccard_distance_matrix.tsv",
        sep="\t",
        index=True,
    )

    coords, eigenvalues = classical_mds(distance_matrix=distance_df.values, n_components=2)

    coords_df = pd.DataFrame(
        {
            "isolate_id": ko_df.index.tolist(),
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
        }
    )
    coords_df.to_csv(
        path_or_buf=out_dir / "ko_pcoa_coordinates.tsv",
        sep="\t",
        index=False,
    )

    variance_explained = (eigenvalues / eigenvalues.sum()) * 100
    variance_df = pd.DataFrame(
        {
            "axis": [f"PCoA{i + 1}" for i in range(len(variance_explained))],
            "eigenvalue": eigenvalues,
            "variance_explained_pct": variance_explained,
        }
    )
    variance_df.to_csv(
        path_or_buf=out_dir / "ko_pcoa_variance_explained.tsv",
        sep="\t",
        index=False,
    )

    save_pcoa_plot(
        coords_df=coords_df,
        variance_df=variance_df,
        output_path=out_dir / "ko_pcoa_plot.png",
    )

    membership_df, summary_df = build_intersection_tables(ko_df=ko_df)
    membership_df.to_csv(
        path_or_buf=out_dir / "ko_isolate_membership_long.tsv",
        sep="\t",
        index=False,
    )
    summary_df.to_csv(
        path_or_buf=out_dir / "ko_intersection_summary.tsv",
        sep="\t",
        index=False,
    )

    save_intersection_plot(
        summary_df=summary_df,
        output_path=out_dir / "ko_top_intersections.png",
        top_n=args.top_n_intersections,
    )


if __name__ == "__main__":
    main()
