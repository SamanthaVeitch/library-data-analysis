#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--input", default="data/METABRIC_RNA_Mutation.csv",
                    help="Kaggle breast cancer dataset")
parser.add_argument("--output", default="figures/histologic_grade_clustering.png")
parser.add_argument("--seed", default=0, type=int)
args = parser.parse_args()

# hypothesis: gene expression patterns can separate breast cancers by
# histologic grade (1, 2, or 3).
# note: gene expression values are normalized z-scores relative to
# reference samples.
# approach: i) count num. cases per histologic grade,
# ii) use MiniBatchKMeans with k=3 to cluster patients based on gene
# expression, iii) compare cluster assignments with actual tumor grades.
# figures: i) bar plot of histologic grades, ii) cluster sizes,
# iii) histologic grade distribution within each cluster.

def main():
    np.random.seed(args.seed)
    df = pd.read_csv(args.input) # read data into pandas DataFrame
    ### organizing columns into 3 categories: i) metadata, ii) gene expression, and
    # iii) gene mutation
    
    mutation_col_suffix = '_mut'
    col_idx_gene_expr_start = -1 # first column of gene expression
    col_idx_mut_start = -1 # first column of gene mutation
    col_idx_histologic_grade = -1 # target
    
    for i,col in enumerate(df.columns):
        if col == 'neoplasm_histologic_grade':
            col_idx_histologic_grade = i
        if col == 'brca1':
            col_idx_gene_expr_start = i
        if col.endswith(mutation_col_suffix):
            col_idx_mut_start = i
            break # exit loop
    # verify col indices were discovered
    assert (
        all([
            col_idx_gene_expr_start != -1,
            col_idx_mut_start != -1,
            col_idx_histologic_grade != -1
        ])
    )

    ### histologic grade distribution
    target = df.iloc[:, col_idx_histologic_grade]
    target_counts = target.value_counts().sort_index()

    print("Histologic grade distribution:")
    print(target_counts)

    ### gene expression data
    x = df.iloc[:, col_idx_gene_expr_start:col_idx_mut_start]
    assert (x.isna().sum().sum() == 0)

    ### MiniBatchKMeans clustering
    model = MiniBatchKMeans(
        n_clusters=3,
        random_state=args.seed
    )

    clusters = model.fit_predict(x)

    ### compare clusters with histologic grade
    results = pd.DataFrame({
        'Cluster': clusters,
        'Histologic Grade': target
    })

    # remove patients without a known histologic grade
    results = results.dropna(subset=['Histologic Grade'])

    # count each grade within each cluster
    cluster_grade_counts = pd.crosstab(
        results['Cluster'],
        results['Histologic Grade']
    )

    # calculate percentage of each grade within each cluster
    cluster_grade_percent = cluster_grade_counts.div(
        cluster_grade_counts.sum(axis=1),
        axis=0
    ) * 100

    ### visualization
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
    ax = ax.flatten()

    # histologic grade distribution
    ax[0].bar(x=target_counts.index,
              height=target_counts.values, color='black')
    ax[0].set_title("Histologic Grade Distribution", loc="left")
    ax[0].set_xlabel("Histologic Grade")
    ax[0].set_ylabel("Count")

    # number of patients in each cluster
    cluster_counts = pd.Series(clusters).value_counts().sort_index()
    ax[1].bar(x=cluster_counts.index,
              height=cluster_counts.values, color='black')
    ax[1].set_title("Cluster Sizes", loc="left")
    ax[1].set_xlabel("Cluster")
    ax[1].set_ylabel("Count")

    # histologic grade distribution within each cluster
    cluster_grade_percent.plot(
        kind='bar',
        stacked=True,
        ax=ax[2]
    )

    ax[2].set_title("Histologic Grade by Cluster", loc="left")
    ax[2].set_xlabel("Cluster")
    ax[2].set_ylabel("Percentage (%)")
    ax[2].legend(title="Histologic Grade")

    for i in range(len(ax)):
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)

    plt.tight_layout()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_path)
    plt.close(fig)


if __name__ == "__main__":
    main()