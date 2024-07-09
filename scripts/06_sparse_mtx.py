import warnings
import numpy as np
import pandas as pd
from ete3 import Tree
import anndata as ad
import scipy as sp
from multiprocessing import Pool
from itertools import repeat
from scipy.sparse import coo_matrix

warnings.simplefilter(action='ignore', category=FutureWarning)

def create_node_ids(subset, dataset, base_path):
    """Create node IDs from a Newick tree file."""
    tree_path = base_path + f"{dataset}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"
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
    data_path = base_path + f'{dataset}/cellstatesappli/RNAmatrix.tsv'
    data = pd.read_csv(data_path, delimiter='\t', index_col=0, nrows=1)
    return data.columns.tolist()

def get_cells_barcodes(node, base_path):
    """Retrieve cell barcodes associated with a node."""
    dataset = node.split('-')[1]
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
    return row_indices, col_indices, values

if __name__ == "__main__":
    base_path = '/home/celia/Documents/VisualStudioCode/'
    dataset = 'WERBA'
    samples = ['P01', 'P02','P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15',
               'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27']
    output_file = '/media/celia/data/cell_mtx_werba.npz'
    num_processes = 3

    nodes = create_dataset_node_list(dataset, base_path)
    cells = create_dataset_cell_list(dataset, base_path)

    cells_set = set(cells)

    with Pool(processes=num_processes) as pool:
        results = pool.starmap(process_node, zip(range(len(nodes)), nodes, repeat(cells_set), repeat(cells), repeat(base_path)))

    row_indices, col_indices, values = [], [], []
    for result in results:
        row_indices.extend(result[0])
        col_indices.extend(result[1])
        values.extend(result[2])

    num_nodes = len(nodes)
    num_cells = len(cells)
    node_cell_matrix = coo_matrix((values, (row_indices, col_indices)), shape=(num_nodes, num_cells))
    sp.sparse.save_npz(output_file, node_cell_matrix)
