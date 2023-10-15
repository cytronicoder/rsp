import dash
from dash import dcc, html
from dash.dependencies import Input, Output

from scripts.tsne import generate_tsne
from scripts.simulation import plot_simulated_cells
from scripts.rsp import generate_polygon

# Initialize the Dash app
app = dash.Dash(__name__)

# Define the app layout
app.layout = html.Div([
    html.H1("RSP Analysis"),
    dcc.Graph(id='rsp-plot'),  # Placeholder for the plot
    html.Label('Marker Gene Percentage'),
    dcc.Slider(
        id='marker-percentage-slider',
        min=1, # can't be 0
        max=100,
        value=50,
        marks={i: f"{i}%" for i in range(0, 101, 10)},
        step=0.1
    )
])

# Define the callback to update the plot based on the slider value
@app.callback(
    Output('rsp-plot', 'figure'),
    Input('marker-percentage-slider', 'value')
)
def update_rsp_plot(percentage_value):
    expression_percentage = percentage_value / 100
    coordinates, is_expressing = plot_simulated_cells(
        num_points=1000, expression_percentage=expression_percentage, distribution="biased", sigma=0.3, seed=42
    )
    # Assuming rsp_analysis returns a figure. If not, modify accordingly.
    fig = generate_polygon(coordinates, is_expressing)
    return fig

# Run the Dash app
if __name__ == '__main__':
    app.run_server(debug=True)
