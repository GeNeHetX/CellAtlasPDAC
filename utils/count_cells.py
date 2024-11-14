from pathlib import Path
import scipy.sparse as sp
from ete3 import Tree

# Dataset information
datasets = {
    'peng': ['T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22'],
    'steele': ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11A", "11B", "12", "13", "15", "16"],
    'lin': ["P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "MET01", "MET02", "MET03", "MET04", "MET05", "MET06"],
    'hwang': ['003_U15', '004_U18', '007_U2', '008_T22', '009_T20', '010T_U12', '011_U9', '2083_T4', '2100_T2', '2223_T17',
              '2229_T18', '2276_U4', '2364_U5', '2376_U1', '2443_U10', '2462_T16', '2490_U11', '2498_U7', '2507_T25', '2523_U3',
              '2540_T6', '2591_U13', '2603_U6', '2626_U16', '2634_T5', '2664_U14', '2667_T23', '2668_T7', '2675_T13', 'MGH2076_T9',
              'MGH2101_T8', 'MGH2381_T19', 'MGH2675_T13', 'MGHR1_T14', 'MGHR2_T15', 'MGHR3_T3', 'MGHR5_T1', 'MGHR6_T10', 'MGHR7_T12', 'MGHR8_T24',
              'MGHR9_T21', 'MGHR11_T11', 'MGHR16_U17', 'MGHR17_U8'],
    'zhang': ['Case1.YF', 'Case1.ZY', 'Case2.YF', 'Case2.ZC', 'Case2.ZY', 'Case3.YF', 'Case3.ZY', 'Case4.ZY'],
    'werba': ['P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15',
               'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27']
}

def load_data_():
    # Helper functions
    def create_node_ids(subset, dataset):
        tree_path = f"/home/celia/Documents/VisualStudioCode/{dataset.upper()}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"
        t = Tree(tree_path, format=1)
        return [f"{node.name}_{subset}-{dataset}" for node in t.traverse() if not node.is_root()]

    def create_dataset_node_list(dataset):
        return [node for subset in datasets[dataset] for node in create_node_ids(subset, dataset)]

    sparse_matrices = {}
    node_lists = {}
    

    sparse_matrix_dir = Path('/media/celia/data/cell_mtx/')
    for dataset in datasets:
        # Load sparse matrix
        sparse_matrix_path = sparse_matrix_dir / f'cell_mtx_{dataset}.npz'
        if sparse_matrix_path.exists():
            sparse_matrices[dataset] = sp.load_npz(sparse_matrix_path)
        else:
            print(f"Warning: Sparse matrix file not found for dataset {dataset} at {sparse_matrix_path}")

        # Load nodes
        node_lists[dataset] = create_dataset_node_list(dataset)

    return sparse_matrices, node_lists

def get_nb_cells(node, sparse_matrices, node_lists):
    dataset = node.split('-')[-1].lower()
    nodes = node_lists[dataset]
    node_index = nodes.index(node)
    matching_indices = sparse_matrices[dataset].row == node_index
    return matching_indices.sum()

sparse_matrices_1, node_lists_1 = load_data_()

get_nb_cells('C1_T2-peng', sparse_matrices_1, node_lists_1)