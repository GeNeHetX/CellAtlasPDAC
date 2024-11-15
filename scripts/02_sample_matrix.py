import pandas as pd
from pathlib import Path
import os

sample_output_dir = Path('/home/celia/Documents/VisualStudioCode/samples/')
dataset = 'steele'

RNA_matrix_path = Path(f"/home/celia/Documents/VisualStudioCode/{dataset.upper()}/cellstatesappli/RNAmatrix.tsv")
RNA_matrix = pd.read_csv(RNA_matrix_path, sep="\t", index_col=0)

sample_separator = '_'
sample_position = 2

# Extract sample IDs from column names
sample_ids = list(set(item.split(sample_separator)[sample_position] for item in RNA_matrix.columns))
print(sample_ids)

# Filter the RNA matrix for each sample ID and save the results
for sample_id in sample_ids:
    columns_to_keep = RNA_matrix.columns[RNA_matrix.columns.str.contains(f'PDAC_TISSUE_{sample_id}_')]
    RNA_filter = RNA_matrix[columns_to_keep]

    # Filter genes with 0 expression in all cells
    RNA_filter = RNA_filter[RNA_filter.sum(axis=1) != 0]

    # Filter out cells with 0 expression in all genes
    RNA_filter = RNA_filter.loc[:, (RNA_filter != 0).any(axis=0)]

    # Save the filtered matrix to a new file
    filtered_path = Path(f"/home/celia/Documents/VisualStudioCode/{dataset.upper()}/cellstatesappli/RNAmatrix_{sample_id}.tsv")
    RNA_filter.to_csv(filtered_path, sep="\t")

# # Create directories for each sample ID
# for sample_id in sample_ids:
#     try:
#         os.makedirs(f"/cellstatesappli/RNAmatrix_{sample_id}", exist_ok=True)
#     except Exception as e:
#         print(f"Error creating directory for {sample_id}: {e}")


# file_path = sample_output_dir / f'{dataset}_samples.txt'
# with open(file_path, 'w') as f:
#     for sample in sample_ids:
#         f.write(f"{sample}\n")

# print(f"Written {len(sample_ids)} samples to {file_path}")