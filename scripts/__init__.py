"""
The `scripts` module provides utility functions for:

1. Saving Plotly figures to various file formats.
2. Plotting t-SNE coordinates with highlights based on expression data.
3. Generating t-SNE 2D coordinates from DGE files and performing clustering using DBSCAN.
4. Simulating cellular data by generating points within a unit circle, and plotting based on gene expression.

Key Functions
----------

- `save_plot(fig: go.Figure, filename: str) -> None`:
    Save the given Plotly figure to a file.

- `plot(tsne_coordinates: array, is_expressing: array, title: str) -> None`:
    Plot t-SNE coordinates highlighting specific points.

- `generate_tsne(dge_file: str, output_file: str, marker_gene: str, target_cluster: int, epsilon: int, minpts: int, dev: bool) -> array`:
    Process DGE files and generate t-SNE coordinates, and perform clustering using DBSCAN.

- `plot_simulated_cells(num_points: int, expression_percentage: float, distribution: str, sigma: float, seed: int, display: bool) -> Tuple[array, array]`:
    Simulate cellular data by generating points and plotting based on gene expression. Two variations of this function are provided.

Examples, detailed parameters, and returns are documented within the individual function docstrings.

"""
