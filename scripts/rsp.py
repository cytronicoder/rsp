"""
This module offers utilities for performing gene analysis using both t-SNE and RSP (Radar Scanning Plot) representations.

The primary functionality allows users to:

1. Generate a polygonal representation of data distributions given certain coordinates.
2. Conduct a comprehensive gene analysis by generating t-SNE plots and RSP plots for a specified marker gene and target cluster.
"""

import numpy as np
import plotly.graph_objects as go

from scripts.tsne import generate_tsne


def generate_polygon(coordinates, is_expressing, theta_bound=[0, 2 * np.pi]):
    """
    Generates a polygonal representation of data distributions using given coordinates.

    Parameters:
    - coordinates (np.ndarray): A 2D array of x and y coordinates.
    - is_expressing (np.ndarray): Boolean array indicating the expressing status of each coordinate.
    - theta_bound (list, optional): Boundaries for the angle theta. Default is [0, 2 * np.pi].

    Returns:
    - rsp_fig (plotly.graph_objects.Figure): A Plotly figure representing the generated polygon.

    Example:
    >>> coords = np.array([[1,2],[3,4],[5,6]])
    >>> expressing = np.array([True, False, True])
    >>> fig = generate_polygon(coords, expressing)
    >>> fig.show()
    """
    resolution = 1000
    theta_start, theta_end = theta_bound
    angle_step = (theta_end - theta_start) / resolution

    is_expressing = is_expressing.flatten()
    foreground_coordinates = coordinates[is_expressing]
    background_coordinates = coordinates[~is_expressing]

    total_expressing_cells = sum(is_expressing)

    differences = []

    # Edge case: the particular gene is expressed in all cells
    #            usually, this applies to some mitocondrial genes
    if len(background_coordinates) == 0:
        differences = [0] * resolution
    else:
        for theta in np.linspace(theta_start, theta_end, resolution):
            # Compute projections
            foreground_projection = (
                np.cos(theta) * foreground_coordinates[:, 0]
                + np.sin(theta) * foreground_coordinates[:, 1]
            )
            background_projection = (
                np.cos(theta) * background_coordinates[:, 0]
                + np.sin(theta) * background_coordinates[:, 1]
            )

            min_val = min(foreground_projection.min(), background_projection.min())
            max_val = max(foreground_projection.max(), background_projection.max())

            norm_denom = max_val - min_val
            foreground_projection = (foreground_projection - min_val) / norm_denom
            background_projection = (background_projection - min_val) / norm_denom

            foreground_hist_values, _ = np.histogram(
                foreground_projection, bins=resolution, range=(0, 1)
            )
            background_hist_values, _ = np.histogram(
                background_projection, bins=resolution, range=(0, 1)
            )

            # Compute CDFs
            foreground_cdf = np.cumsum(foreground_hist_values)
            foreground_cdf = foreground_cdf / foreground_cdf[-1]

            background_cdf = np.cumsum(background_hist_values)
            background_cdf = background_cdf / background_cdf[-1]

            difference = foreground_cdf - background_cdf
            differences.append(np.sum(np.abs(difference)))

    # Construct the polygon
    angles = []
    radii = []

    for theta, diff in zip(
        np.linspace(theta_start, theta_end, resolution), differences
    ):
        h = diff / (np.sin(angle_step / 2) * resolution)
        angles.extend([theta, theta + angle_step])
        radii.extend([h, h])

    # Normalize radii
    radii = [r / total_expressing_cells for r in radii]

    segment_areas = []
    for i in range(0, len(radii) - 1, 2):
        theta_diff = angles[i + 1] - angles[i]
        segment_area = 0.5 * theta_diff * (radii[i] ** 2 + radii[i + 1] ** 2)
        segment_areas.append(segment_area)

    polygon_area = np.sum(segment_areas)

    radii_np = np.array(radii)
    biggest_radius = np.max(radii_np)

    circle_area = np.pi * biggest_radius**2
    area_difference = np.abs(polygon_area - circle_area)
    roundness = area_difference / circle_area

    # Plotly polar plot for the polygon
    rsp_fig = go.Figure()
    rsp_fig.add_trace(
        go.Scatterpolar(
            r=radii,
            theta=np.degrees(angles),
            fill="toself",
            name="Polygon",
            line=dict(color="blue"),
            opacity=0.5,
        )
    )

    rsp_fig.update_layout(
        title="Polygon Representation",
        polar=dict(
            radialaxis=dict(visible=True, range=[0, 1]),
            angularaxis=dict(
                direction="counterclockwise",
                showticklabels=True,
                ticks="outside",
                rotation=0,
            ),
        ),
        annotations=[
            go.layout.Annotation(
                x=0,
                y=0,
                text=f"Area: {polygon_area:.2f}",
                showarrow=False,
                font=dict(size=16),
            )
        ],
        showlegend=True,
    )

    return rsp_fig, polygon_area, roundness


def gene_analysis(
    dge_file,
    marker_gene=None,
    target_cluster=None,
    theta_bound=[0, 2 * np.pi],
    debug=False,
):
    """
    Performs gene analysis by generating t-SNE and RSP plots for a given marker gene and target cluster.

    Parameters:
    - dge_file (str): Path to the DGE (Differential Gene Expression) file.
    - marker_gene (str, optional): The marker gene of interest. Default is None.
    - target_cluster (str, optional): The target cluster to focus on. Default is None.
    - theta_bound (list, optional): Boundaries for the angle theta used in the RSP plot. Default is [0, 2 * np.pi].

    Returns:
    - tsne_fig (plotly.graph_objects.Figure): A Plotly figure representing the t-SNE plot.
    - rsp_fig (plotly.graph_objects.Figure): A Plotly figure representing the RSP plot.

    Example:
    >>> tsne_fig, rsp_fig = gene_analysis('path/to/dge/file.txt', marker_gene='GeneA', target_cluster='Cluster1')
    >>> tsne_fig.show()
    >>> rsp_fig.show()
    """
    tsne_coordinates, is_expressing, tsne_fig = generate_tsne(
        dge_file, marker_gene=marker_gene, target_cluster=target_cluster, debug=debug
    )

    rsp_fig, rsp_area, roundness = generate_polygon(
        tsne_coordinates, is_expressing, theta_bound=theta_bound
    )

    rsp_fig.update_layout(
        title=f"RSP plot with marker gene '{marker_gene}'",
    )

    return tsne_fig, rsp_fig, rsp_area, roundness
