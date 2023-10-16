import numpy as np
import plotly.graph_objects as go
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import webbrowser
from threading import Timer

from scripts.tsne import generate_tsne
from scripts.simulation import plot_simulated_cells

marker_gene = input("Please enter the name of the marker gene: ")

coordinates, is_expressing = generate_tsne(
    dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
    marker_gene=marker_gene,
)

# coordinates, is_expressing = plot_simulated_cells(
#     num_points=1000, expression_percentage=0.80, distribution="biased", sigma=0.3
# )

initial_theta = np.radians(float(input("Please enter an initial angle (in degrees): ")))

app = dash.Dash(__name__)


def compute_z_scores(coordinates, theta):
    x, y = coordinates[:, 0], coordinates[:, 1]
    z_scores = x * np.cos(theta) + y * np.sin(theta)
    return z_scores - np.min(z_scores)


# App layout
app.layout = html.Div(
    [
        dcc.Graph(id="histogram-graph"),
        dcc.Graph(id="difference-graph"),
        dcc.Slider(
            id="theta-slider",
            min=initial_theta - np.pi / 2,
            max=initial_theta + np.pi / 2,
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


@app.callback(
    [
        Output("histogram-graph", "figure"),
        Output("difference-graph", "figure"),
        Output("theta-display", "children"),
    ],
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

    fg_hist_values, fg_hist_bins = np.histogram(fg_z_scores)
    bg_hist_values, bg_hist_bins = np.histogram(bg_z_scores)
    diff_values = fg_hist_values - bg_hist_values

    diff_fig = go.Figure()
    diff_fig.add_trace(
        go.Bar(
            x=fg_hist_bins[:-1], y=diff_values, name="Difference", marker_color="green"
        )
    )
    diff_fig.update_layout(
        title=f"Difference Histogram at θ = {np.degrees(theta):.2f}°"
    )

    return fig, diff_fig, f"θ = {np.degrees(theta):.2f}°"


if __name__ == "__main__":
    Timer(1, lambda: webbrowser.open("http://127.0.0.1:8050/")).start()
    app.run_server(debug=True, use_reloader=False)
