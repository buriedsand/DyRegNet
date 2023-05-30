import sys
import pandas as pd

ugrn_path = sys.argv[1]
H3K27ac_RP_path = sys.argv[2]
TR_RP_path = sys.argv[3]
output_file = sys.argv[4]

# Load data
ugrn = pd.read_csv(ugrn_path, names=["tr", "gene"])
TR_RP = pd.read_csv(TR_RP_path, index_col=0).sort_index()
H3K27ac_RP = pd.read_csv(H3K27ac_RP_path, index_col=0, names=["H3K27ac"]).sort_index()

# Normalize data
TR_RP_normalized = (TR_RP - TR_RP.min()) / (TR_RP.max() - TR_RP.min())
TR_RP_normalized.fillna(0, inplace=True)
H3K27ac_RP_normalized = (H3K27ac_RP - H3K27ac_RP.min()) / (
    H3K27ac_RP.max() - H3K27ac_RP.min()
)

# Repeat the H3K27ac_RP_normalized to match the number of columns in TR_RP_normalized
H3K27ac_RP_normalized_repeated = pd.concat(
    [H3K27ac_RP_normalized] * TR_RP_normalized.shape[1], axis=1
)
H3K27ac_RP_normalized_repeated.columns = TR_RP_normalized.columns

# Perform element-wise multiplication
result = TR_RP_normalized.multiply(H3K27ac_RP_normalized_repeated)

# Prepare data for grouping
result.reset_index(inplace=True)
result["index"] = result["index"].map(lambda x: x.split(":")[1])

# Group by 'index' column and select the maximum value for each group
grouped_result = result.groupby("index").max()

# Create a series with multiindex 'gene' and 'tr'
weights_series = grouped_result.stack()
weights_series.index.names = ["gene", "tr"]

# Set the index of ugrn to be a MultiIndex with 'gene' and 'tr'
ugrn = ugrn.set_index(["gene", "tr"])

# Update the 'weight' column using the values from the weights_series
ugrn["weight"] = 0  # Initialize the 'weight' column
ugrn["weight"].update(weights_series)

# Reset the index
ugrn = ugrn.reset_index()

# Rearrange columns
ugrn = ugrn[["tr", "gene", "weight"]]

# ugrn = ugrn[ugrn["weight"] > 0]
# print(ugrn["weight"].max())
print(ugrn)

ugrn.to_csv(output_file, index=False)
