"""
TSNE generation script and related functions.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.metrics import pairwise_distances

dge_data = None
expression_matrix = None
tsne_coordinates = None
cluster_labels = None


def plot_k_distance_graph(k):
    """
    Generate k-distance graph for t-SNE coordinates; useful in determining
    the value of epsilon for DBSCAN.

    Parameters:
    - k (int): The value of k for computing the k-distance graph.
    """
    global tsne_coordinates

    dists = pairwise_distances(tsne_coordinates)
    k_dists = np.sort(dists, axis=1)[:, k]
    sorted_k_dists = np.sort(k_dists)

    plt.plot(sorted_k_dists)
    plt.title("k-distance graph")
    plt.xlabel(f"Points sorted with {k}th nearest distances")
    plt.ylabel(f"{k}th nearest distances")
    plt.show()


def plot(dge_file, output_file, epsilon=4, minpts=40, dev=False, marker_gene=None):
    """
    Generate t-SNE 2D coordinates from DGE file and cluster using DBSCAN.

    Parameters:
    - dge_file (str): The path to the DGE file.
    - output_file (str): The path to the output file.
    - epsilon (int): The epsilon value for DBSCAN.
    - minpts (int): The minpts value for DBSCAN.
    - dev (bool): Whether to print debug information.
    - marker_gene (str): The name of the marker gene to highlight.
    """
    global dge_data, expression_matrix, tsne_coordinates, cluster_labels

    if os.path.isfile(output_file):
        # Read the coordinates from the output file if it exists
        tsne_coordinates = pd.read_csv(output_file).values
    else:
        # Read the DGE file and generate t-SNE coordinates
        dge_data = pd.read_csv(dge_file, sep=None, engine="python")
        expression_matrix = dge_data.iloc[:, 1:].values.T.astype(float)

        if dev:
            print(
                f"Loaded {expression_matrix.shape[1]} genes in {expression_matrix.shape[0]} cells."
            )

        tsne = TSNE()
        tsne_coordinates = tsne.fit_transform(expression_matrix)
        pd.DataFrame(tsne_coordinates).to_csv(
            output_file, index=False, header=["X", "Y"]
        )

    if dev:
        plot_k_distance_graph(minpts)

    # Calculate cluster IDs using DBSCAN
    dbscan = DBSCAN(eps=epsilon, min_samples=minpts)
    cluster_labels = dbscan.fit_predict(tsne_coordinates)
    cluster_labels[cluster_labels != -1] = cluster_labels[cluster_labels != -1] + 1

    if dev:
        # Print number of noise points and number of points in each cluster
        print(
            "Number of noise points: {}".format(
                len(tsne_coordinates[cluster_labels == -1])
            )
        )
        for i in range(1, max(cluster_labels) + 1):
            print(
                "Number of points in cluster {}: {}".format(
                    i, len(tsne_coordinates[cluster_labels == i])
                )
            )

    if marker_gene:
        # TODO: Highlight a specific gene in red (foreground) and all other genes in grey (background)
        pass
    else:
        unique_labels = np.unique(cluster_labels)
        for label in unique_labels:
            mask = cluster_labels == label
            plt.scatter(
                tsne_coordinates[mask, 0],
                tsne_coordinates[mask, 1],
                s=1,
                alpha=0.1,
                label=f"Cluster {label}" if label != -1 else "Noise",
            )
            if label != -1:
                cluster_center = tsne_coordinates[mask].mean(axis=0)
                plt.annotate(
                    str(label),
                    cluster_center,
                    fontsize=10,
                    ha="center",
                    va="center",
                    color="black",
                    weight="bold",
                )

    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.title("t-SNE plot with DBSCAN clustering")
    plt.legend(loc="best", markerscale=5)
    plt.show()
