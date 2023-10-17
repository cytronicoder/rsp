import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import atexit
import os
import tempfile
import base64
import webbrowser
from threading import Timer
import shutil
from scripts.rsp import gene_analysis

# Initialize the Dash app
app = dash.Dash(__name__)

styles = {
    "input": {
        "width": "20%",
        "padding": "10px",
        "margin": "5px",
        "display": "inline-block",
        "border": "1px solid #ccc",
        "borderRadius": "4px",
        "boxSizing": "border-box",
    },
    "button": {
        "padding": "10px 15px",
        "margin": "10px 0",
        "border": "none",
        "borderRadius": "4px",
        "cursor": "pointer",
    },
    "upload": {
        "border": "1px solid #ccc",
        "display": "inline-block",
        "padding": "6px 12px",
        "cursor": "pointer",
        "margin": "5px 0",
        "borderRadius": "4px",
    },
}

# Define the app layout
app.layout = html.Div(
    [
        html.H1("RSP Analysis"),
        html.Div(
            [
                dcc.Upload(
                    id="upload-dge-file",
                    children=html.Button("Upload DGE File", style=styles["button"]),
                    style=styles["upload"],
                ),
                dcc.Input(
                    id="gene-input",
                    type="text",
                    placeholder="Enter Gene Name",
                    style=styles["input"],
                ),
                dcc.Input(
                    id="cluster-input",
                    type="number",
                    placeholder="Enter Cluster",
                    style=styles["input"],
                ),
                html.Button(
                    "Generate Plots", id="generate-button", style=styles["button"]
                ),
                html.Div(
                    [
                        dcc.Graph(id="tsne-plot", style={"width": "50%"}),
                        dcc.Graph(id="rsp-plot", style={"width": "50%"}),
                    ],
                    style={"display": "flex", "justifyContent": "space-between"},
                ),
            ],
        ),
    ]
)

temp_dir = tempfile.mkdtemp()  # Create a dedicated temporary directory for the app


@app.callback(
    Output("upload-dge-file", "children"),
    [Input("upload-dge-file", "contents")],
    prevent_initial_call=True,
)
def update_upload_button_text(file_contents):
    if not file_contents:
        return html.Button("Upload DGE File", style=styles["button"])

    dge_file_path = save_dge_to_temp_file(file_contents, "dge")
    return f"File uploaded to {dge_file_path}. Click again to change."


@app.callback(
    [Output("tsne-plot", "figure"), Output("rsp-plot", "figure")],
    [Input("generate-button", "n_clicks")],  # Listening to the button click event
    [
        # These are now 'State' not 'Input' because we only want to read their values, not trigger the callback
        dash.dependencies.State("upload-dge-file", "contents"),
        dash.dependencies.State("gene-input", "value"),
        dash.dependencies.State("cluster-input", "value"),
    ],
)
def update_plots(n_clicks, dge_file_content, gene_name, cluster):
    if not n_clicks or dge_file_content is None:
        return dash.no_update, dash.no_update

    dge_file_path = save_dge_to_temp_file(dge_file_content, "dge")

    tsne_fig, rsp_fig, _ = gene_analysis(
        dge_file_path, marker_gene=gene_name, target_cluster=cluster, debug=False
    )

    return tsne_fig, rsp_fig


def save_dge_to_temp_file(content, filename):
    content_type, content_string = content.split(",")
    decoded = base64.b64decode(content_string)
    dge_file_path = os.path.join(temp_dir, f"{filename}_temp.txt")

    with open(dge_file_path, "wb") as f:
        f.write(decoded)

    # Log the saved path
    print(f"Saved DGE file to {dge_file_path}")

    return dge_file_path


# Cleanup function
def cleanup_temp_files():
    try:
        print(f"Deleting temp directory: {temp_dir}")
        shutil.rmtree(temp_dir)  # Remove the entire directory
    except Exception as e:
        print(f"Error deleting temp directory {temp_dir}: {e}")


atexit.register(cleanup_temp_files)

# Run the Dash app
if __name__ == "__main__":
    Timer(1, lambda: webbrowser.open("http://127.0.0.1:8050/")).start()
    app.run_server(debug=True, use_reloader=False)
