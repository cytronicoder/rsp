import os
from scripts.rsp import gene_analysis
from scripts.util import get_genes, save_plot

os.makedirs("plots", exist_ok=True)

genes = get_genes(dge_file="data/GSM2906447_NeonatalHeart_dge.txt", target_cluster=1)

for gene in genes:
    print(f"Generating plots for {gene}...")
    tsne_fig, rsp_fig = gene_analysis(
        dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
        marker_gene=gene,
        target_cluster=1,
    )

    for fig in [tsne_fig, rsp_fig]:
        os.makedirs(f"plots/{gene}", exist_ok=True)
        save_plot(fig, f"plots/{gene}/{fig.layout.title.text}.png")
