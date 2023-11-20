import pandas as pd
from scipy.stats import pearsonr
import itertools

# Load your data: replace 'your_data.csv' with your file path
data = pd.read_csv("data/MCA1_lung.csv")
data.set_index("GENE", inplace=True)

# Potential genes to compare
potential_genes = [
    # ...
]

# Calculate Pearson Correlation Coefficients
correlations = {}
p_values = {}

# Iterate over all possible pairs of genes
for gene1, gene2 in itertools.combinations(potential_genes, 2):
    if gene1 in data.index and gene2 in data.index:
        corr, p_value = pearsonr(data.loc[gene1], data.loc[gene2])
        correlations[(gene1, gene2)] = corr
        p_values[(gene1, gene2)] = p_value

# Convert to DataFrame for easier handling
correlation_df = pd.DataFrame.from_dict(
    correlations, orient="index", columns=["Pearson_Correlation"]
)
p_value_df = pd.DataFrame.from_dict(p_values, orient="index", columns=["P_Value"])

# Combine both DataFrames
result_df = pd.concat([correlation_df, p_value_df], axis=1)

# Filtering for significant results
# Adjust the threshold as necessary
threshold = 0.0
significant_results = result_df[
    (result_df["Pearson_Correlation"].abs() >= threshold)
    & (result_df["P_Value"] < 0.05)
]

# Save results to a CSV file
significant_results.to_csv("significant_correlations.csv")

print("Analysis complete. Results saved to 'significant_correlations.csv'")
