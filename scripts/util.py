import os
import pandas as pd

from scripts.tsne import generate_tsne


def get_genes(dge_file, target_cluster=None):
    """
    Get genes from the DGE file. If target_cluster is specified, use DBSCAN to
    isolate and return only the genes present in that cluster after t-SNE reduction.

    Parameters:
    - dge_file (str): The path to the DGE file.
    - target_cluster (int, optional): The target cluster to isolate. Defaults to None.

    Returns:
    - list: List of genes.
    """

    # Use the generate_tsne function to get t-SNE coordinates and clustering information
    _, is_expressing, _ = generate_tsne(dge_file, target_cluster=target_cluster)

    # Read the DGE file to get gene names
    split_filename = os.path.splitext(dge_file)[0]
    if os.path.isfile(f"{split_filename}.dge.parquet"):
        dge_data = pd.read_parquet(f"{split_filename}.dge.parquet")
    else:
        dge_data = pd.read_csv(dge_file, sep=None, engine="python")

    gene_names = dge_data.index.tolist()

    print(is_expressing)

    # If a target cluster is specified, return only genes from that cluster
    if target_cluster is not None and is_expressing is not None:
        genes_in_target_cluster = [
            gene for i, gene in enumerate(gene_names) if is_expressing[i]
        ]
        return genes_in_target_cluster

    return gene_names


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
