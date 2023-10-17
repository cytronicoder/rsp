import os
import pandas as pd

from scripts.rsp import gene_analysis
from scripts.util import get_genes, get_gene_info, save_plot

# os.makedirs("plots", exist_ok=True)


# Define a function to collect rows
def collect_rows(rows, new_row):
    rows.append(new_row)
    return rows


# List to store rows
rows = []

genes = get_genes(dge_file="data/GSM2906447_NeonatalHeart_dge.txt", target_cluster=1)

for gene in genes:
    print(f"Reading {gene}...")
    # Get the gene information
    info = get_gene_info(
        dge_file="data/GSM2906447_NeonatalHeart_dge.txt", target_gene=gene
    )

    # Get the RSP Area from the gene_analysis function
    _, _, rsp_area = gene_analysis(
        dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
        marker_gene=gene,
        target_cluster=1,
    )

    # Collect data row by row
    new_row = {
        "Gene Name": info[0],
        "Coverage (%)": info[1],
        "Mean Expression": info[2],
        "Total Expression": info[3],
        "RSP Area": rsp_area,
    }
    rows = collect_rows(rows, new_row)

    # print(f"Generating plots for {gene}...")
    # tsne_fig, rsp_fig = gene_analysis(
    #     dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
    #     marker_gene=gene,
    #     target_cluster=1,
    # )

    # for fig in [tsne_fig, rsp_fig]:
    #     os.makedirs(f"plots/{gene}", exist_ok=True)
    #     save_plot(fig, f"plots/{gene}/{fig.layout.title.text}.png")


# Convert the list of rows into a DataFrame
master_df = pd.DataFrame(rows)

# Save the master DataFrame to a CSV file
master_df.to_csv("master_file.csv", index=False)
