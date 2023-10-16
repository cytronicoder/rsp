import os
import numpy as np
import pandas as pd

from scripts.tsne import generate_tsne


def get_genes(dge_file, target_cluster=None):
    """
    Get genes from the DGE file. If target_cluster is specified, filter and
    return only the genes that are expressed in that cluster.

    Parameters:
    - dge_file (str): The path to the DGE file.
    - target_cluster (int, optional): The target cluster to filter by. Defaults to None.

    Returns:
    - list: List of genes.
    """

    # Use the generate_tsne function to get t-SNE coordinates and clustering information
    _, is_expressing_cells, _ = generate_tsne(dge_file, target_cluster=target_cluster)

    # Read the DGE file to get gene names and expression matrix
    split_filename = os.path.splitext(dge_file)[0]
    if os.path.isfile(f"{split_filename}.dge.parquet"):
        dge_data = pd.read_parquet(f"{split_filename}.dge.parquet")
    else:
        dge_data = pd.read_csv(dge_file, sep=None, engine="python")

    gene_names = dge_data.index.values
    expression_matrix = dge_data.values

    # Filter the expression matrix to only include cells from the target cluster
    if is_expressing_cells is not None:
        reduced_matrix = expression_matrix[:, is_expressing_cells]
    else:
        reduced_matrix = expression_matrix

    # Get the list of genes that are expressed in any of the cells in the target cluster
    expressed_genes = [
        gene
        for gene, expression in zip(gene_names, reduced_matrix)
        if np.any(expression > 0)
    ]

    return expressed_genes


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
