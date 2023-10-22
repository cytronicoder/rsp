import pandas as pd
import scanpy as sc
from sklearn.cluster import DBSCAN


def generate_tsne(
    dge_file,
    output_file=None,
    marker_gene=None,
    target_cluster=None,
    epsilon=4,
    minpts=40,
    debug=False,
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
    - debug (bool): Whether to print debug information.

    Returns:
    - filtered_tsne_coordinates (numpy.ndarray): A filtered set of t-SNE coordinates. If a marker gene or target cluster is specified, this contains only the coordinates for cells expressing the gene or belonging to the target cluster, respectively.
    - is_expressing (numpy.ndarray or None): A boolean array indicating which cells (of the filtered set) are expressing the specified marker gene. If no marker gene is specified, this is None.
    - fig (plotly.graph_objects.Figure): A Plotly figure object visualizing the t-SNE plot with either DBSCAN clusters or highlighted cells based on marker gene expression.
    """

    if debug:
        print(f"Running in debug mode!")

    # Load the data
    adata = sc.read_csv(dge_file, delimiter="\t").T

    # Perform PCA before t-SNE
    sc.pp.pca(adata, n_comps=30)

    # Run t-SNE on the principal components
    sc.tl.tsne(adata, use_rep="X_pca")

    # DBSCAN clustering
    dbscan = DBSCAN(eps=epsilon, min_samples=minpts)
    cluster_labels = dbscan.fit_predict(adata.obsm["X_tsne"])
    cluster_labels[cluster_labels != -1] = cluster_labels[cluster_labels != -1] + 1

    # Add cluster labels to adata
    adata.obs["DBSCAN"] = cluster_labels

    # Save the t-SNE coordinates to the output file if specified
    if output_file:
        pd.DataFrame(adata.obsm["X_tsne"], columns=["X", "Y"]).to_csv(
            output_file, index=False
        )

    # Prepare for plotting
    if marker_gene:
        if marker_gene in adata.var_names:
            if target_cluster is not None:
                adata = adata[adata.obs['DBSCAN'] == target_cluster].copy()
            
            gene_expression = adata[:, marker_gene].X
            is_expressing = gene_expression > 0
        else:
            print(f"Marker gene '{marker_gene}' not found in the DGE data.")
            return None, None, None
    elif target_cluster is not None:
        adata = adata[adata.obs["DBSCAN"] == target_cluster].copy()
        is_expressing = None
    else:
        is_expressing = None

    # Plot the data using scanpy's built-in functions
    if marker_gene:
        sc.pl.tsne(
            adata,
            color=marker_gene,
            show=True,
            title=f"t-SNE plot with {marker_gene} highlighted",
        )
    else:
        sc.pl.tsne(
            adata, color="DBSCAN", show=False, title="t-SNE plot with DBSCAN clustering"
        )

    return (
        adata.obsm["X_tsne"],
        is_expressing,
        None,
    )  # The third None is because scanpy's plot doesn't return a fig by default
