import warnings
import numpy as np
import pandas as pd
from ete3 import Tree
import anndata as ad
import scipy as sp
from multiprocessing import Pool
from itertools import repeat
from scipy.sparse import coo_matrix
from tqdm import tqdm
import gc
from pathlib import Path
import os

warnings.simplefilter(action='ignore', category=FutureWarning)

def create_node_ids(subset, dataset, base_path):
    """Create node IDs from a Newick tree file."""
    tree_path = base_path + f"{dataset.upper()}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"
    t = Tree(tree_path, format=1)
    return [f"{node.name}_{subset}-{dataset}" for node in t.traverse() if not node.is_root()]

def create_dataset_node_list(dataset, base_path):
    """Create a list of node IDs for a dataset."""
    dataset_node_list = []
    for subset in samples:
        dataset_node_list.extend(create_node_ids(subset, dataset, base_path))
    return dataset_node_list

def create_dataset_cell_list(dataset, base_path):
    """Create a list of cell barcodes for a dataset."""
    data_path = base_path + f'{dataset.upper()}/cellstatesappli/RNAmatrix.tsv'
    data = pd.read_csv(data_path, delimiter='\t', index_col=0, nrows=1)
    return data.columns.tolist()

def get_cells_barcodes(node, base_path):
    """Retrieve cell barcodes associated with a node."""
    dataset = node.split('-')[1].upper()
    node_parts = node.split('-')[0]
    subset = '_'.join(node_parts.split('_')[1:])
    node_name = node_parts.split('_')[0]

    tree_path = base_path + f"{dataset}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"
    t = Tree(tree_path, format=1)
    node = t.search_nodes(name=node_name)
    leaves = [leaf.name for leaf in node[0].iter_leaves()]

    data_path = base_path + f'{dataset}/cellstatesappli/RNAmatrix_{subset}.tsv'
    data = pd.read_csv(data_path, delimiter='\t', index_col=0).T

    cluster_path = base_path + f'{dataset}/cellstatesappli/RNAmatrix_{subset}/optimized_clusters.txt'
    clust = np.loadtxt(cluster_path, dtype=np.int64)

    adata = ad.AnnData(X=data.values)
    adata.var_names = data.columns
    adata.obs_names = data.index
    adata.obs['cluster'] = clust

    clusters = [int(i[1:]) for i in leaves]
    cell_list = []
    for clst in clusters:
        bdata = adata[adata.obs['cluster'] == clst]
        cell_list.extend(bdata.obs_names.tolist())
    return cell_list


def process_node(node_index, node, cells_set, cells, base_path):
    """Process a single node to find its cell indices and values."""
    row_indices = []
    col_indices = []
    values = []
    for cell in set(get_cells_barcodes(node, base_path)).intersection(cells_set):
        row_indices.append(node_index)
        col_indices.append(cells.index(cell))
        values.append(1)
    return np.array(row_indices), np.array(col_indices), np.array(values)


if __name__ == "__main__":
    base_path = '/home/celia/Documents/VisualStudioCode/'
    sample_dir = '/home/celia/Documents/VisualStudioCode/samples/'
    dataset = 'peng'
    output_file = f'/media/celia/data/cell_mtx/cell_mtx_{dataset}_data.npz'
    num_processes = 2

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"Output file already exists for sample: {dataset}. Skipping processing.")
        exit() 

    sample_dir = Path(sample_dir)
    file_path = sample_dir / f'{dataset}_samples.txt'
    
    if file_path.exists():
        with open(file_path, 'r') as f:
            samples = [line.strip() for line in f.readlines()]
    else:
        print(f"Warning: File {file_path} not found!")

    nodes = create_dataset_node_list(dataset, base_path)
    cells = create_dataset_cell_list(dataset, base_path)

    cells_set = set(cells)

    print(f"Starting multiprocessing with {len(nodes)} nodes and {len(cells)} cells...")

    # Initialize the progress bar
    pbar = tqdm(total=len(nodes), desc="Processing Nodes", ncols=100)

    def update_progress_bar(result):
        """Callback to collect results and update the progress bar."""
        row_indices.append(result[0])  # Append row indices from each process
        col_indices.append(result[1])  # Append col indices
        values.append(result[2])  # Append values
        pbar.update()  # Update the progress bar

    # Preallocate result containers as lists of numpy arrays
    row_indices, col_indices, values = [], [], []

    # Use multiprocessing
    with Pool(processes=num_processes) as pool:
        for node_index, node in enumerate(nodes):
            pool.apply_async(process_node, (node_index, node, cells_set, cells, base_path), callback=update_progress_bar)
        
        pool.close()
        pool.join()

    pbar.close()
    print("Multiprocessing completed.")

    # Concatenate all results to avoid memory issues with large lists
    row_indices = np.concatenate(row_indices)
    col_indices = np.concatenate(col_indices)
    values = np.concatenate(values)

    num_nodes = len(nodes)
    num_cells = len(cells)

    # Construct sparse matrix using the COO format
    node_cell_matrix = coo_matrix((values, (row_indices, col_indices)), shape=(num_nodes, num_cells))
    
    np.savez_compressed(
        output_file,
        sparse_matrix=node_cell_matrix,
        cell_index_map=np.array(list(cells), dtype=object),
        node_index_map=np.array(list(nodes), dtype=object)
    )


    # Free up memory by manually invoking garbage collection
    del row_indices, col_indices, values
    gc.collect()

    print(f"Matrix saved to {output_file}")