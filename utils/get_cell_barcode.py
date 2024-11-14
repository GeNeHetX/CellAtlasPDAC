import numpy as np
import os
from pathlib import Path
from scipy.sparse import csr_matrix

datasets_name = ['peng' , 'steele', 'lin' , 'hwang', 'zhang', 'werba']
cell_mtx_dir = Path('/media/celia/data/cell_mtx/')

node_lists = {}
cells_lists = {}
sparse_matrices = {}

for dataset in datasets_name:
    sparse_data_path = cell_mtx_dir / f'cell_mtx_{dataset}_data.npz'
    if sparse_data_path.exists():
        data = np.load(sparse_data_path, allow_pickle=True)
        sparse_matrices[dataset] = csr_matrix(data['sparse_matrix'].item())
        node_lists[dataset] = data['node_index_map']
        cells_lists[dataset] = data['cell_index_map']
    else:
        print(f"Warning: Sparse matrix file not found for dataset {dataset} at {sparse_data_path}")

def get_cell_barcodes(node):
    dataset = node.split('-')[-1]
    node_index_map = node_lists[dataset]
    sparse_mtx = sparse_matrices[dataset]
    cell_index_map = cells_lists[dataset]
    if node not in node_index_map:
        raise ValueError(f"Node '{node}' not found in the index map")

    row_idx = np.where(node_index_map == node)[0][0]
    cells_indices = sparse_mtx[row_idx].nonzero()[1]
    cell_barcodes = cell_index_map[cells_indices]
    return cell_barcodes

print(get_cell_barcodes('C1_T2-peng'))

# node = 'I1_P05-lin'

# file_path = f'/media/celia/data/cell_mtx/cell_mtx_lin_data.npz'

# # Check if file exists
# if not os.path.exists(file_path):
#     raise FileNotFoundError(f"Data file not found: {file_path}")

# # Load the data file
# data = np.load(file_path, allow_pickle=True)

# # Extract data from the npz file
# try:
#     sparse_matrix = data['sparse_matrix'].item()
#     cell_index_map = data['cell_index_map']
#     node_index_map = data['node_index_map']
# except KeyError as e:
#     raise KeyError(f"Missing key in data file: {e}")

# sparse_matrix = csr_matrix(sparse_matrix)

# # Check if the node exists in the index map
# if node not in node_index_map:
#     raise ValueError(f"Node '{node}' not found in the index map")

# # Get the index of the node
# row_idx = np.where(node_index_map == node)[0][0]

# # # Extract indices of cells
# cells_indices = sparse_matrix[row_idx].nonzero()[1]
# cell_barcodes = cell_index_map[cells_indices]


