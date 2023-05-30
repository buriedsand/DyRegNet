import sys
import pandas as pd
from pathlib import Path

input_files = sys.argv[1:-1]
output_file = sys.argv[-1]


# Initialize an empty DataFrame
result = pd.DataFrame()

# Iterate over the CSV files
for file in input_files:
    factor_name = Path(file).stem

    # Read the CSV file into a DataFrame
    df = pd.read_csv(file, header=None, names=["Gene", factor_name])

    # Merge the DataFrames, using the "Gene" column as the key
    if result.empty:
        result = df
    else:
        result = result.merge(df, on="Gene", how="outer")

# Sort the DataFrame by the "Gene" column
result = result.sort_values(by="Gene")

# Reset the index and drop the old index
result = result.reset_index(drop=True)

# Set the index as the "Gene" column
result = result.set_index("Gene")

# Remove index name
result = result.rename_axis(None)

# Sort columns alphabetically
result = result.sort_index(axis=1)

# Save the aggregated DataFrame to a new CSV file
result.to_csv(output_file)
