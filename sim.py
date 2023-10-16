import dash
from dash import dcc, html
from dash.dependencies import Input, Output

from scripts.tsne import generate_tsne
from scripts.simulation import plot_simulated_cells
from scripts.rsp import generate_polygon

# Initialize the Dash app
app = dash.Dash(__name__)

# Define the app layout
app.layout = html.Div(
    [
        html.H1("RSP Analysis"),
        html.Div(
            [
                html.Div(
                    [
                        dcc.Graph(id="tsne-plot", style={"width": "50%"}),
                        dcc.Graph(id="rsp-plot", style={"width": "50%"}),
                    ],
                    style={"display": "flex", "justifyContent": "space-between"},
                ),
                html.Div(
                    [
                        html.Label("Marker Gene Percentage"),
                        dcc.Slider(
                            id="marker-percentage-slider",
                            min=1,
                            max=100,
                            value=50,
                            marks={i: f"{i}%" for i in range(0, 101, 10)},
                            step=0.1,
                        ),
                    ],
                ),
                html.Div(
                    [
                        html.Label("Optional Seed"),
                        dcc.Input(
                            id="seed-input",
                            type="number",
                            placeholder="Enter a seed (optional)",
                            style={"margin-left": "10px"},
                        ),
                    ],
                    style={"display": "flex"},
                ),
                html.Div(
                    [
                        html.Label("Distribution Type"),
                        dcc.RadioItems(
                            id="distribution-type",
                            options=[
                                {"label": "Biased", "value": "biased"},
                                {"label": "Even", "value": "even"},
                            ],
                            value="biased",
                            style={"margin-left": "10px"},
                        ),
                    ],
                    style={"display": "flex"},
                ),
            ],
            style={"width": "100%", "display": "flex", "flex-direction": "column", "gap": "10px"},
        )
    ]
)


# Define the callback to update the plots based on the slider, seed value, and distribution type
@app.callback(
    [Output("tsne-plot", "figure"), Output("rsp-plot", "figure")],
    [
        Input("marker-percentage-slider", "value"),
        Input("seed-input", "value"),
        Input("distribution-type", "value"),
    ],
)
def update_plots(percentage_value, seed_value, distribution_type):
    expression_percentage = percentage_value / 100
    seed = seed_value if seed_value is not None else 42
    distribution_type = distribution_type if distribution_type is not None else "biased"
    coordinates, is_expressing, fig1 = plot_simulated_cells(
        num_points=1000,
        expression_percentage=expression_percentage,
        distribution=distribution_type,
        sigma=0.3,
        seed=seed,
    )
    fig2 = generate_polygon(coordinates, is_expressing)
    return fig1, fig2


# Run the Dash app
if __name__ == "__main__":
    app.run_server(debug=True)
