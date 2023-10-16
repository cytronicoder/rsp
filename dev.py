from scripts.tsne import generate_tsne
from scripts.simulation import plot_simulated_cells
from scripts.rsp import generate_polygon, gene_analysis
from scripts.util import save_plot

tsne_fig, rsp_fig = gene_analysis(
    dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
    marker_gene="Actc1",
    target_cluster=1,
)

for fig in [tsne_fig, rsp_fig]:
    fig.show()

# coordinates, is_expressing = generate_tsne(
#     dge_file="data/GSM2906447_NeonatalHeart_dge.txt",
#     marker_gene="Actc1",
#     target_cluster=1,
# )

# coordinates, is_expressing, _ = plot_simulated_cells(
#     num_points=1000,
#     expression_percentage=0.50,
#     distribution="biased",
#     sigma=0.3,
#     display=False,
# )

# fig = generate_polygon(coordinates, is_expressing)
# save_plot(fig, "rsp.png")

# from 0 to 1 in intervals of 0.1
# for i in range(1, 10, 1):
#     expression_percentage = i / 10
#     coordinates, is_expressing, _ = plot_simulated_cells(
#         num_points=1000,
#         expression_percentage=expression_percentage,
#         distribution="biased",
#         seed=12,
#         display=True,
#     )

#     fig = generate_polygon(coordinates, is_expressing)
#     fig.show()

#     print(input(f"{expression_percentage * 100}%; press enter to continue..."))


# # range from 10% to 100% in intervals of 10%
# for i in range(11, 100, 11):
#     expression_percentage = i / 100
#     coordinates, is_expressing, _ = plot_simulated_cells(
#         num_points=1000,
#         expression_percentage=expression_percentage,
#         distribution="biased",
#         sigma=0.3,
#         # seed=42,
#         display=False,
#     )

#     fig = generate_polygon(coordinates, is_expressing)
#     fig.show()

#     print(input("Press enter to continue..."))
