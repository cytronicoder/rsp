"""
TSNE generation script and related functions.
"""

import os
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.metrics import pairwise_distances
import plotly.graph_objects as go
import plotly.express as px

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

    fig = px.line(
        x=np.arange(len(sorted_k_dists)),
        y=sorted_k_dists,
        labels={"x": "Points", "y": f"{k}th nearest distances"},
    )
    fig.update_layout(title="k-distance graph")
    fig.show()


def generate_tsne(
    dge_file,
    output_file=None,
    marker_gene=None,
    target_cluster=None,
    epsilon=4,
    minpts=40,
    dev=False,
):
    """
    Generate t-SNE 2D coordinates from DGE file and cluster using DBSCAN.

    Parameters:
    - dge_file (str): The path to the DGE file.
    - output_file (str): The path to the output file; defaults to
      the same directory as the DGE file with a .tsne.csv extension.
    - marker_gene (str): The name of the marker gene to highlight.
    - target_cluster (int): The target cluster to highlight.
    - epsilon (int): The epsilon value for DBSCAN.
    - minpts (int): The minpts value for DBSCAN.
    - dev (bool): Whether to print debug information.
    """
    global dge_data, expression_matrix, tsne_coordinates, cluster_labels
    print(f"Running in dev mode: {dev}")

    split_filename = os.path.splitext(dge_file)[0]

    if output_file is None:
        output_file = split_filename + ".tsne.csv"
        print(f"Defaulting output file directory to {output_file}.")

    # Read the DGE file
    if os.path.isfile(f"{split_filename}.dge.parquet"):
        # Reading cache
        dge_data = pd.read_parquet(f"{split_filename}.dge.parquet")
    else:
        dge_data = pd.read_csv(dge_file, sep=None, engine="python")

        # Caching
        dge_data.to_parquet(f"{split_filename}.dge.parquet")

    if os.path.isfile(output_file):
        # Read the coordinates from the output file if it exists
        tsne_coordinates = pd.read_csv(output_file).values
    else:
        # Generate t-SNE coordinates
        expression_matrix = dge_data.values.T.astype(float)

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
        print(input("Press Enter to continue..."))

    # Calculate cluster IDs using DBSCAN
    dbscan = DBSCAN(eps=epsilon, min_samples=minpts)
    cluster_labels = dbscan.fit_predict(tsne_coordinates)
    cluster_labels[cluster_labels != -1] = cluster_labels[cluster_labels != -1] + 1

    fig = go.Figure()

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

    # Filter for returning later
    filtered_tsne_coordinates = tsne_coordinates

    if marker_gene:
        if marker_gene in dge_data.index:
            gene_expression = dge_data.loc[marker_gene].values
            is_expressing = gene_expression > 0
            foreground_mask = is_expressing

            if target_cluster is not None:
                cluster_mask = cluster_labels == target_cluster
                is_expressing = is_expressing & cluster_mask
                filtered_tsne_coordinates = tsne_coordinates[cluster_mask]
                is_expressing = is_expressing[cluster_mask]
                foreground_mask = foreground_mask & cluster_mask
                background_mask = cluster_mask & ~foreground_mask
            else:
                background_mask = ~foreground_mask

            fig.add_trace(
                go.Scatter(
                    x=tsne_coordinates[background_mask, 0],
                    y=tsne_coordinates[background_mask, 1],
                    mode="markers",
                    marker=dict(color="lightgray", size=5, opacity=0.75),
                    name="Background",
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=tsne_coordinates[foreground_mask, 0],
                    y=tsne_coordinates[foreground_mask, 1],
                    mode="markers",
                    marker=dict(color="red", size=5),
                    name=f"Foreground ({marker_gene})",
                )
            )
        else:
            print(f"Marker gene '{marker_gene}' not found in the DGE data.")
            return None, None
    elif target_cluster is not None:
        is_expressing = cluster_labels == target_cluster
        filtered_tsne_coordinates = tsne_coordinates[is_expressing]
        mask = cluster_labels == target_cluster

        fig.add_trace(
            go.Scatter(
                x=tsne_coordinates[mask, 0],
                y=tsne_coordinates[mask, 1],
                mode="markers",
                marker=dict(size=5),
                name=f"Cluster {target_cluster}",
            )
        )
    else:
        # Plot all clusters if no marker gene is specified
        is_expressing = None  # No marker gene, so no foreground/background
        unique_labels = np.unique(cluster_labels)
        for label in unique_labels:
            mask = cluster_labels == label
            fig.add_trace(
                go.Scatter(
                    x=tsne_coordinates[mask, 0],
                    y=tsne_coordinates[mask, 1],
                    mode="markers",
                    marker=dict(size=5, opacity=0.1 if label == -1 else 1.0),
                    name=f"Cluster {label}" if label != -1 else "Noise",
                )
            )

    fig.update_layout(
        title="t-SNE plot with DBSCAN clustering"
        if not marker_gene
        else f"t-SNE plot with {marker_gene} highlighted"
    )
    fig.show()

    return filtered_tsne_coordinates, is_expressing


def plot(tsne_coordinates, is_expressing, title="t-SNE Plot"):
    """
    Plot t-SNE coordinates highlighting points based on the is_expressing array.

    Parameters:
    - tsne_coordinates (array): 2D array of t-SNE coordinates.
    - is_expressing (array): Boolean array indicating which points to highlight.
    - title (str): The title for the plot.
    """

    if tsne_coordinates is None or is_expressing is None:
        print("No t-SNE coordinates or is_expressing array provided.")
        return

    # Points that are expressing
    foreground_mask = np.array(is_expressing)
    background_mask = ~foreground_mask

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=tsne_coordinates[background_mask, 0],
            y=tsne_coordinates[background_mask, 1],
            mode="markers",
            marker=dict(color="lightgray", size=5),
            name="Background",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=tsne_coordinates[foreground_mask, 0],
            y=tsne_coordinates[foreground_mask, 1],
            mode="markers",
            marker=dict(color="red", size=5),
            name="Foreground",
        )
    )

    fig.update_layout(title=title)
    fig.show()
