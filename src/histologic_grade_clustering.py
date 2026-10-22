#!/usr/bin/env python3
"""
Analyze the relationship between gene expression and histologic tumor grade.

This script uses MiniBatchKMeans to group breast cancer samples from the
METABRIC dataset based on normalized gene-expression measurements.

Hypothesis:
    Gene-expression patterns can distinguish breast cancer tumors
    associated with different histologic grades (1, 2, and 3).

Methods:
    1. Load the METABRIC dataset and identify gene-expression columns.
    2. Cluster samples using MiniBatchKMeans with K=3.
    3. Compare cluster assignments with known histologic grades.
    4. Generate plots showing grade counts, cluster sizes, and the
       histologic grade distribution within each cluster.

The histologic grades are not used to fit the clustering model.
"""

import argparse
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score


def main():
    """Cluster METABRIC samples and visualize histologic grade patterns."""

    # Parse command-line arguments.
    parser = argparse.ArgumentParser(
        description="Cluster METABRIC gene expression by histologic grade."
    )
    parser.add_argument(
        "--input",
        default="data/METABRIC_RNA_Mutation.csv",
        help="Path to the METABRIC dataset"
    )
    parser.add_argument(
        "--output",
        default="figures/histologic_grade_clustering.png",
        help="Path to save the resulting figure"
    )
    parser.add_argument(
        "--seed",
        default=0,
        type=int,
        help="Random seed for reproducible clustering"
    )
    args = parser.parse_args()

    # Load the METABRIC dataset into a pandas DataFrame.
    df = pd.read_csv(args.input, low_memory=False)

    # The dataset contains three column groups:
    # 1. Patient and tumor metadata
    # 2. Normalized gene-expression measurements
    # 3. Gene-mutation information
    #
    # Identify the beginning of the expression and mutation columns,
    # along with the histologic grade column.
    mutation_col_suffix = "_mut"
    col_idx_gene_expr_start = -1
    col_idx_mut_start = -1
    col_idx_histologic_grade = -1

    for i, col in enumerate(df.columns):
        if col == "neoplasm_histologic_grade":
            col_idx_histologic_grade = i

        if col == "brca1":
            col_idx_gene_expr_start = i

        if col.endswith(mutation_col_suffix):
            col_idx_mut_start = i
            break

    # Verify that all required columns were identified.
    if not all([
        col_idx_gene_expr_start != -1,
        col_idx_mut_start != -1,
        col_idx_histologic_grade != -1
    ]):
        raise ValueError("Required METABRIC columns were not found.")

    # Extract known histologic grades (1, 2, or 3).
    # Missing grades are excluded automatically from value_counts().
    target = df.iloc[:, col_idx_histologic_grade]
    target_counts = target.value_counts().sort_index()

    print("Histologic grade distribution:")
    print(target_counts)

    # Extract the normalized gene-expression features.
    # Each row represents a tumor sample and each column a gene.
    # Mutation columns and clinical metadata are excluded.
    x = df.iloc[:, col_idx_gene_expr_start:col_idx_mut_start]

    # MiniBatchKMeans requires complete numerical input data.
    if x.isna().any().any():
        raise ValueError("Gene-expression data contains missing values.")

    # Group samples into three clusters using gene expression only.
    # The random seed makes the clustering reproducible.
    model = MiniBatchKMeans(
        n_clusters=3,
        random_state=args.seed
    )

    clusters = model.fit_predict(x)

    # Associate each sample's cluster with its known histologic grade.
    results = pd.DataFrame({
        "Cluster": clusters,
        "Histologic Grade": target
    })

    # Samples with unknown grades are still included in clustering,
    # but cannot be included in the grade comparison.
    results = results.dropna(subset=["Histologic Grade"])

    # Count the number of tumors of each grade in each cluster.
    cluster_grade_counts = pd.crosstab(
        results["Cluster"],
        results["Histologic Grade"]
    )

    # Convert grade counts into percentages within each cluster.
    # Each row sums to 100%.
    cluster_grade_percent = cluster_grade_counts.div(
        cluster_grade_counts.sum(axis=1),
        axis=0
    ) * 100

    print("\nHistologic grade percentages by cluster:")
    print(cluster_grade_percent.round(2))

    # Count all samples assigned to each cluster, including those
    # without a known histologic grade.
    cluster_counts = pd.Series(clusters).value_counts().sort_index()

    print("\nCluster sizes:")
    print(cluster_counts)

    # Calculate the total variance explained by the PCA transformation.
    pca = PCA(n_components=2)
    x_pca = pca.fit_transform(x)
    explained_variance = pca.explained_variance_ratio_.sum() * 100
    print("\nExplained Variance:")
    print(f"{explained_variance:.2f}%")

    # Calculate the silhouette score to evaluate clustering quality.
    score = silhouette_score(x, clusters)
    print("\nSilhouette Score:")
    print(score)

    # Create three plots to summarize the clustering analysis.
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
    ax = ax.flatten()

    # Plot 1: Number of tumors in each histologic grade.
    ax[0].bar(
        x=target_counts.index,
        height=target_counts.values,
        color="black"
    )
    ax[0].set_title("Histologic Grade Distribution", loc="left")
    ax[0].set_xlabel("Histologic Grade")
    ax[0].set_ylabel("Count")

    # Plot 2: Scatter plot of tumors by cluster assignment.
    # Plot tumors colored by their MiniBatchKMeans cluster.
    ax[1].scatter = ax[1].scatter(
        x_pca[:, 0],
        x_pca[:, 1],
        c=clusters,
        cmap="viridis",
        s=12,
        alpha=0.6
    )

    ax[1].set_xlabel(
        f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}%"
        f" variance explained)"
    )

    ax[1].set_ylabel(
        f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}%"
        f" variance explained)"
    )
    ax[1].set_title("MiniBatchKMeans Breast Cancer Clusters")

    ax[1].legend(
        handles=ax[1].scatter.legend_elements()[0],
        labels=["Cluster 0", "Cluster 1", "Cluster 2"]
    )

    # Plot 3: Percentage of each histologic grade within each cluster.
    cluster_grade_percent.plot(
        kind="bar",
        stacked=True,
        ax=ax[2]
    )
    ax[2].set_title("Histologic Grade by Cluster", loc="left")
    ax[2].set_xlabel("Cluster")
    ax[2].set_ylabel("Percentage (%)")
    ax[2].legend(title="Histologic Grade")

    # Remove the top and right plot borders for readability.
    for axis in ax:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    # Create the output directory if it does not already exist.
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save the final figure and close it to release resources.
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)


if __name__ == "__main__":
    main()
