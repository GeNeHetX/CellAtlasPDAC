import os
import numpy as np
from scipy.sparse import load_npz

gene_mtx_folder_path = '/media/celia/data/gene_mtx'
node_name = 'I25_004_U18-hwang' 

# Extract dataset, sample, and node from the node_name
dataset = node_name.split('-')[-1]
print(f"Dataset: {dataset}")
sample = '_'.join(node_name.split('-')[0].split('_')[1:])
print(f"Sample: {sample}")
node = node_name.split('-')[0].split('_')[0]
print(f"Node: {node}")

# Construct paths
mtx_sample_path = f"{gene_mtx_folder_path}/{dataset}/{sample}"
file_path = f'{mtx_sample_path}/{dataset}_{sample}_data.npz'

# Check if file exists
if not os.path.exists(file_path):
    raise FileNotFoundError(f"Data file not found: {file_path}")

# Load the data file
data = np.load(file_path, allow_pickle=True)

# Extract data from the npz file
try:
    sparse_matrix_over = data['sparse_matrix_over'].item()
    sparse_matrix_under = data['sparse_matrix_under'].item()
    gene_index_map_over = dict(data['gene_index_map_over'])
    gene_index_map_under = dict(data['gene_index_map_under'])
    node_index_map = dict(data['node_index_map'])
except KeyError as e:
    raise KeyError(f"Missing key in data file: {e}")

# Check if the node exists in the index map
if node not in node_index_map:
    raise ValueError(f"Node '{node}' not found in the index map")

# Get the index of the node
row_idx = node_index_map[node]

# Extract indices of over-expressed genes
over_gene_indices = sparse_matrix_over[row_idx].nonzero()[1]
over_expressed_genes = np.array(list(gene_index_map_over.keys()))[over_gene_indices]
print(f"Over-expressed genes for {node_name}: {over_expressed_genes}")

# Extract indices of under-expressed genes
under_gene_indices = sparse_matrix_under[row_idx].nonzero()[1]
under_expressed_genes = np.array(list(gene_index_map_under.keys()))[under_gene_indices]
print(f"Under-expressed genes for {node_name}: {under_expressed_genes}")


