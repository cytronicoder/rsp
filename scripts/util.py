import os
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from scripts.tsne import generate_tsne


def get_genes(dge_file, target_cluster=None):
    """
    Get genes from the DGE file. If target_cluster is specified, use DBSCAN to
    isolate and return only the genes present in that cluster.

    Parameters:
    - dge_file (str): The path to the DGE file.
    - target_cluster (int, optional): The target cluster to isolate. Defaults to None.

    Returns:
    - list: List of genes.
    """

    # Read the DGE file
    split_filename = os.path.splitext(dge_file)[0]

    if os.path.isfile(f"{split_filename}.dge.parquet"):
        dge_data = pd.read_parquet(f"{split_filename}.dge.parquet")
    else:
        dge_data = pd.read_csv(dge_file, sep=None, engine="python")
        dge_data.to_parquet(f"{split_filename}.dge.parquet")

    if target_cluster is not None:
        tsne_coordinates, _, _ = generate_tsne(dge_file, target_cluster=target_cluster)
        dbscan = DBSCAN()
        cluster_labels = dbscan.fit_predict(tsne_coordinates)
        target_cell_indices = np.where(cluster_labels == target_cluster)[0].tolist()

        filtered_data = dge_data[dge_data.columns[target_cell_indices].tolist()]
        expressed_genes_mask = filtered_data.sum(axis=1) > 1
        expressed_genes = filtered_data.index[expressed_genes_mask].tolist()

        return expressed_genes

    # If no target cluster specified, return all genes
    return dge_data.index.tolist()


def save_plot(fig, filename):
    """
    Save the given Plotly figure to a file.

    Parameters:
    - fig (go.Figure): The Plotly figure to save.
    - filename (str): The name of the file to save the figure to.
      Supported extensions: .html, .png, .jpg or .jpeg

    Example:
    >>> save_plot(fig, 'my_plot.html')
    """
    file_ext = filename.split(".")[-1]

    if file_ext == "html":
        fig.write_html(filename)
    elif file_ext in ["png", "jpg", "jpeg"]:
        fig.write_image(filename)
    else:
        print(f"Unsupported file format: {file_ext}. Use one of: html, png, jpg, jpeg.")
