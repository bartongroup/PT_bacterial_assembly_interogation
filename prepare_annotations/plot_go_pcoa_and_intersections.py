#!/usr/bin/env python3

"""
Generate GO-based PCoA outputs and GO intersection summaries.
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
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate GO PCoA outputs and GO intersection summaries."
    )
    parser.add_argument(
        "--go_matrix_tsv",
        required=True,
        help="Path to GO presence/absence matrix TSV.",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for GO PCoA and intersection outputs.",
    )
    parser.add_argument(
        "--top_n_intersections",
        type=int,
        default=20,
        help="Number of top GO intersections to plot.",
    )
    return parser.parse_args()


def load_matrix(matrix_tsv: Path) -> pd.DataFrame:
    """Load a binary matrix indexed by isolate_id."""
    df = pd.read_csv(matrix_tsv, sep="\t", dtype=str, low_memory=False)
    if "isolate_id" not in df.columns:
        raise ValueError("Expected 'isolate_id' column in matrix.")
    df = df.set_index("isolate_id")
    df = df.apply(pd.to_numeric)
    return df


def compute_jaccard_distance_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Compute Jaccard distance matrix between isolates."""
    distances = pdist(X=df.values, metric="jaccard")
    distance_matrix = squareform(X=distances)
    return pd.DataFrame(distance_matrix, index=df.index, columns=df.index)


def classical_mds(
    distance_matrix: np.ndarray,
    n_components: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    """Perform classical multidimensional scaling."""
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


def build_intersection_tables(
    df: pd.DataFrame,
    feature_name: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Build feature membership and intersection summary tables."""
    membership_records = []

    for feature_id in df.columns:
        present_isolates = df.index[df[feature_id] > 0].tolist()
        membership_key = "|".join(sorted(present_isolates))

        membership_records.append(
            {
                feature_name: feature_id,
                "n_isolates_present": len(present_isolates),
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
        .agg(n_features=(feature_name, "count"))
        .sort_values(
            by=["n_features", "n_isolates_present", "isolates_present"],
            ascending=[False, False, True],
        )
    )

    return membership_df, summary_df


def save_pcoa_plot(
    coords_df: pd.DataFrame,
    variance_df: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    """Save a PCoA scatter plot."""
    pc1_var = variance_df.loc[0, "variance_explained_pct"]
    pc2_var = variance_df.loc[1, "variance_explained_pct"]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(coords_df["PC1"], coords_df["PC2"])

    for _, row in coords_df.iterrows():
        ax.text(row["PC1"], row["PC2"], row["isolate_id"], fontsize=8)

    ax.set_xlabel(f"PCoA1 ({pc1_var:.2f}%)")
    ax.set_ylabel(f"PCoA2 ({pc2_var:.2f}%)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def save_intersection_plot(
    summary_df: pd.DataFrame,
    output_path: Path,
    top_n: int,
    title: str,
) -> None:
    """Save a bar plot of top intersections."""
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
    ax.bar(range(len(top_df)), top_df["n_features"])
    ax.set_xticks(range(len(top_df)))
    ax.set_xticklabels(top_df["plot_label"], rotation=90)
    ax.set_ylabel("Number of GO terms")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    """Run GO PCoA and intersection workflow."""
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    go_df = load_matrix(Path(args.go_matrix_tsv))

    distance_df = compute_jaccard_distance_matrix(go_df)
    distance_df.to_csv(
        out_dir / "go_jaccard_distance_matrix.tsv",
        sep="\t",
        index=True,
    )

    coords, eigenvalues = classical_mds(distance_df.values, n_components=2)

    coords_df = pd.DataFrame(
        {
            "isolate_id": go_df.index.tolist(),
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
        }
    )
    coords_df.to_csv(out_dir / "go_pcoa_coordinates.tsv", sep="\t", index=False)

    variance_explained = (eigenvalues / eigenvalues.sum()) * 100
    variance_df = pd.DataFrame(
        {
            "axis": [f"PCoA{i + 1}" for i in range(len(variance_explained))],
            "eigenvalue": eigenvalues,
            "variance_explained_pct": variance_explained,
        }
    )
    variance_df.to_csv(
        out_dir / "go_pcoa_variance_explained.tsv",
        sep="\t",
        index=False,
    )

    save_pcoa_plot(
        coords_df=coords_df,
        variance_df=variance_df,
        output_path=out_dir / "go_pcoa_plot.png",
        title="PCoA of isolates by GO term presence/absence",
    )

    membership_df, summary_df = build_intersection_tables(go_df, feature_name="go_term")
    membership_df.to_csv(
        out_dir / "go_isolate_membership_long.tsv",
        sep="\t",
        index=False,
    )
    summary_df.to_csv(
        out_dir / "go_intersection_summary.tsv",
        sep="\t",
        index=False,
    )

    save_intersection_plot(
        summary_df=summary_df,
        output_path=out_dir / "go_top_intersections.png",
        top_n=args.top_n_intersections,
        title="Top GO term isolate intersections",
    )


if __name__ == "__main__":
    main()
