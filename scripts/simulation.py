import numpy as np
import plotly.graph_objects as go


def plot_simulated_cells(
    num_points,
    expression_percentage,
    distribution="even",
    sigma=0.3,
    seed=None,
    display=False,
):
    """
    Generate points with coordinates between -1 and 1 using uniform sampling.
    The point is accepted if it satisfies the condition x^2 + y^2 ≤ 1.
    A percentage of the points are randomly chosen to represent the cells
    expressing the gene, and the points are plotted.

    Parameters:
    - num_points (int): Number of points to generate.
    - expression_percentage (float): Percentage of points representing the cells expressing the gene.
    - distribution (str, optional): Can be 'even' for random distribution or 'biased' for clustered foreground.
                                    Defaults to 'even'.
    - sigma (float, optional): Standard deviation of the Gaussian distribution; defaults to 0.3.
    - seed (int, optional): Random seed for reproducibility.
    - display (bool, optional): Whether to display the plot.

    Returns:
    - coordinates (numpy array): Array containing the generated points.
    - is_expressing (numpy array): Boolean array indicating whether each point represents a cell expressing the gene.
    """

    # Setting the seed for reproducibility
    if seed is not None:
        np.random.seed(seed)

    coordinates = np.zeros((num_points, 2))
    is_expressing = np.zeros(num_points, dtype=bool)

    count = 0
    while count < num_points:
        x, y = np.random.uniform(-1, 1, 2)
        if x**2 + y**2 <= 1:
            coordinates[count] = [x, y]
            count += 1

    # Randomly choose a percentage of points to represent the cells expressing the gene
    expressing_indices = np.random.choice(
        num_points, int(num_points * expression_percentage), replace=False
    )

    if distribution == "biased":
        # Define a random center for the Gaussian distribution
        center_x, center_y = np.random.uniform(-0.5, 0.5, 2)

        for idx in expressing_indices:
            x = center_x + np.random.normal(0, sigma)
            y = center_y + np.random.normal(0, sigma)

            # If the point lies outside the main circle, bring it back within the circle
            while x**2 + y**2 > 1:
                x = center_x + np.random.normal(0, sigma)
                y = center_y + np.random.normal(0, sigma)

            coordinates[idx] = [x, y]

    is_expressing[expressing_indices] = True

    if display:
        # Plotting with Plotly
        fig = go.Figure()

        # Plot for background cells
        fig.add_trace(
            go.Scatter(
                x=coordinates[~is_expressing, 0],
                y=coordinates[~is_expressing, 1],
                mode="markers",
                marker=dict(color="gray", size=5, opacity=0.25),
                name="Background",
            )
        )

        # Plot for expressing cells
        fig.add_trace(
            go.Scatter(
                x=coordinates[is_expressing, 0],
                y=coordinates[is_expressing, 1],
                mode="markers",
                marker=dict(color="red", size=5),
                name="Expressing Cells",
            )
        )

        fig.update_layout(
            title=f"{distribution.capitalize()} Distribution with {expression_percentage*100}% Expressing Cells",
            xaxis_title="X",
            yaxis_title="Y",
            width=600,
            height=600,
        )

        fig.show()

    return coordinates, is_expressing
