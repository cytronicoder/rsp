import numpy as np
import matplotlib.pyplot as plt


def compute_z_scores(coordinates, theta):
    """
    Compute z scores for given coordinates and angle theta.

    Parameters:
    - coordinates (np.array): Nx2 array of x, y coordinates.
    - theta (float): angle in radians.

    Returns:
    - z_scores (np.array): z scores for given coordinates and theta.
    """
    x, y = coordinates[:, 0], coordinates[:, 1]
    z_scores = x * np.cos(theta) + y * np.sin(theta)
    return z_scores - np.min(z_scores)


def plot_overlapping_histograms(coordinates, is_expressing_array, theta, bins=30):
    """
    Plot overlapping histograms for foreground and background z-scores.

    Parameters:
    - coordinates (np.array): Nx2 array of x, y coordinates.
    - is_expressing_array (np.array): Boolean array indicating whether each coordinate is in the foreground.
    - theta (float): angle in radians.
    - bins (int): number of bins for the histogram.
    """
    z_scores = compute_z_scores(coordinates, theta)

    fg_z_scores = z_scores[is_expressing_array]
    bg_z_scores = z_scores[~is_expressing_array]

    plt.hist(bg_z_scores, bins=bins, alpha=0.5, label="Background", color="blue")
    plt.hist(fg_z_scores, bins=bins, alpha=0.5, label="Foreground", color="red")
    plt.legend(loc="upper right")
    plt.title(f"Overlapping Histograms at θ = {np.degrees(theta):.2f}°")
    plt.xlabel("Z Score")
    plt.ylabel("Frequency")
    plt.show()
