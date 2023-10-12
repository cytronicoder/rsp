import numpy as np
import plotly.graph_objects as go
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import webbrowser
from threading import Timer

from scripts.tsne import generate_tsne
from scripts.simulation import plot_simulated_cells

# coordinates, is_expressing = generate_tsne(
#     dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
#     marker_gene="Actc1",
#     target_cluster=1,
#     dev=True,
# )

coordinates, is_expressing = plot_simulated_cells(
    num_points=1000, expression_percentage=0.80, distribution="biased", sigma=0.3
)

# Get initial theta from the user
initial_theta = np.radians(float(input("Please enter an initial angle (in degrees): ")))

# Initialize the app
app = dash.Dash(__name__)


# Function to compute z-scores
def compute_z_scores(coordinates, theta):
    x, y = coordinates[:, 0], coordinates[:, 1]
    z_scores = x * np.cos(theta) + y * np.sin(theta)
    return z_scores - np.min(z_scores)


# App layout
app.layout = html.Div(
    [
        dcc.Graph(id="histogram-graph"),
        dcc.Slider(
            id="theta-slider",
            min=initial_theta - np.pi / 2,  # theta - 90 degrees
            max=initial_theta + np.pi / 2,  # theta + 90 degrees
            step=0.01,
            value=initial_theta,
            marks={
                i: f"{int(i*180/np.pi)}°"
                for i in np.linspace(
                    initial_theta - np.pi / 2, initial_theta + np.pi / 2, 9
                )
            },
            tooltip={"placement": "bottom", "always_visible": False},
            updatemode="drag",
        ),
        html.Div(id="theta-display"),
    ]
)


# App callback
@app.callback(
    [Output("histogram-graph", "figure"), Output("theta-display", "children")],
    [Input("theta-slider", "value")],
)
def update_histogram(theta):
    z_scores = compute_z_scores(coordinates, theta)
    fg_z_scores = z_scores[is_expressing]
    bg_z_scores = z_scores[~is_expressing]

    fig = go.Figure()
    fig.add_trace(
        go.Histogram(x=bg_z_scores, name="Background", marker_color="blue", opacity=0.5)
    )
    fig.add_trace(
        go.Histogram(x=fg_z_scores, name="Foreground", marker_color="red", opacity=0.5)
    )
    fig.update_layout(
        barmode="overlay", title=f"Histogram at θ = {np.degrees(theta):.2f}°"
    )

    return fig, f"θ = {np.degrees(theta):.2f}°"


if __name__ == "__main__":
    # Start the server in a different thread to allow for browser open command to execute immediately after
    Timer(1, lambda: webbrowser.open("http://127.0.0.1:8050/")).start()
    app.run_server(debug=True, use_reloader=False)
