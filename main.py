import numpy as np
import matplotlib.pyplot as plt

from scripts.tsne import generate_tsne, plot
from scripts.simulation import plot_simulated_cells

def expression_from_coordinates(coords):
    """
    Convert 2D coordinates to polar coordinates: angles (in radians) and radii.
    """
    radii = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    angles = np.arctan2(coords[:, 1], coords[:, 0]) % (2 * np.pi)  # Ensure angles are [0, 2π]
    return angles, radii

def plot_cumulative_histogram(coordinates, is_expressing, num_bins=100):
    # Convert coordinates to angles and radii
    angles, radii = expression_from_coordinates(coordinates)

    # Split into foreground and background
    foreground_radii = radii[is_expressing]
    background_radii = radii[~is_expressing]

    # Compute the histograms
    fg_hist, fg_edges = np.histogram(foreground_radii, bins=num_bins, density=True)
    bg_hist, bg_edges = np.histogram(background_radii, bins=num_bins, density=True)

    # Calculate the cumulative distributions
    fg_cumulative = np.cumsum(fg_hist) / np.sum(fg_hist)
    bg_cumulative = np.cumsum(bg_hist) / np.sum(bg_hist)

    # Bar width
    bar_width = (fg_edges[1] - fg_edges[0]) * 0.4

    # Create a matplotlib figure
    plt.figure(figsize=(10, 6))

    # Plot bars
    plt.bar(fg_edges[:-1] - bar_width/2, fg_cumulative, width=bar_width, label="Foreground", color="red", align="center")
    plt.bar(bg_edges[:-1] + bar_width/2, bg_cumulative, width=bar_width, label="Background", color="lightgray", alpha=0.75, align="center")

    plt.title("Cumulative Histogram of Radii")
    plt.xlabel("Angle (Radians)")
    plt.ylabel("Cumulative Frequency")
    plt.legend()
    plt.grid(True, which="both", ls="--", c="0.7")
    plt.show()

def plot_difference_histogram(coordinates, is_expressing, num_bins=100):
    # Convert coordinates to angles and radii
    angles, radii = expression_from_coordinates(coordinates)

    # Split into foreground and background
    foreground_radii = radii[is_expressing]
    background_radii = radii[~is_expressing]

    # Compute the histograms
    fg_hist, _ = np.histogram(foreground_radii, bins=num_bins, density=True)
    bg_hist, _ = np.histogram(background_radii, bins=num_bins, density=True)

    # Calculate the difference
    diff_hist = np.abs(fg_hist - bg_hist)

    # Calculate the area of each segment and sum up
    angular_width = 2 * np.pi / num_bins
    segment_areas = 0.5 * angular_width * diff_hist**2
    total_area = np.sum(segment_areas)

    print(f"Total area of the polygon: {total_area:.4f}")

    # Plot
    plt.figure(figsize=(8, 8))
    plt.subplot(111, projection="polar")
    plt.fill(angles, diff_hist, alpha=0.6)
    plt.title("Polygon Representation of Foreground-Background Difference")

    # Show
    plt.show()

    return total_area

coordinates, is_expressing = generate_tsne(
    dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
    marker_gene="Actc1",
    target_cluster=1,
)

plot_cumulative_histogram(coordinates, is_expressing)
plot_difference_histogram(coordinates, is_expressing)
