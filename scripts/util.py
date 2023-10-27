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

    gene_names = dge_data["GENE"].values
    expression_matrix = dge_data.drop(columns=["GENE"]).values

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


def get_gene_info(dge_file, target_gene):
    """
    Get gene information including name, coverage, mean expression, and total expression.

    Parameters:
    - dge_file (str): The path to the DGE file.
    - target_gene (str): The target gene to get information for.

    Returns:
    - tuple: (gene_name, coverage, mean_expression, total_expression) or None if the gene is not found.
    """

    # Load the DGE data
    split_filename = os.path.splitext(dge_file)[0]
    if os.path.isfile(f"{split_filename}.dge.parquet"):
        dge_data = pd.read_parquet(f"{split_filename}.dge.parquet")
    else:
        dge_data = pd.read_csv(dge_file, sep=None, engine="python")

    # Check if the target gene exists in the data
    if target_gene in dge_data["GENE"].values:
        gene_expression_values = (
            dge_data[dge_data["GENE"] == target_gene]
            .drop(columns=["GENE"])
            .values.flatten()
        )

        # Calculate metrics
        gene_name = target_gene

        # Coverage: Percentage of non-zero expressions
        foreground = (gene_expression_values > 0).sum()
        total_samples = len(gene_expression_values)
        coverage = (foreground / total_samples) * 100

        mean_expression = gene_expression_values.mean()
        total_expression = gene_expression_values.sum()

        return (gene_name, coverage, mean_expression, total_expression)

    else:
        print(f"Gene '{target_gene}' not found in the DGE data.")
        return None


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


def migrate_to_new_dge_format(input_path, output_path):
    """
    Process the MCA1.txt file based on specified instructions:
    1. Add a "GENE" column header.
    2. Remove any quotation marks.
    3. Change the delimiter from a space character to a tab.

    Args:
    - input_path: Path to the MCA1.txt file.
    - output_path: Path where the processed file will be saved.

    Returns:
    - Path to the processed file.
    """

    # Reading the input file
    with open(input_path, "r") as file_mca1:
        mca1_content = file_mca1.read()

    # Removing quotation marks
    mca1_content = mca1_content.replace('"', "")

    # Changing the delimiter from space to tab
    mca1_content = mca1_content.replace(" ", "\t")

    # Adding the "GENE" column header and removing "NeonatalHeart_1." prefix
    mca1_lines = mca1_content.splitlines()
    mca1_lines[0] = "GENE\t" + mca1_lines[0].replace("NeonatalHeart_1.", "")
    mca1_processed_content = "\n".join(mca1_lines)

    # Saving the processed content to the output path
    with open(output_path, "w") as output_file:
        output_file.write(mca1_processed_content)

    return output_path
