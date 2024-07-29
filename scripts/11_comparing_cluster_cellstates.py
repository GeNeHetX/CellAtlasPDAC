import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from ete3 import Tree
import scipy.sparse as sp
from tqdm import tqdm

''' CELL BARCODES'''
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
    'werba': ['P01', 'P02', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15',
               'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27']
}

def create_node_ids(subset, dataset):
    tree_path = f"/home/celia/Documents/VisualStudioCode/{dataset.upper()}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"
    t = Tree(tree_path, format=1)
    return [f"{node.name}_{subset}-{dataset}" for node in t.traverse() if not node.is_root()]

def create_dataset_node_list(dataset):
    dataset_node_list = []
    for subset in datasets[dataset]:
        dataset_node_list.extend(create_node_ids(subset, dataset))
    return dataset_node_list

def create_dataset_cell_list(dataset):
    data_path = Path(f'/home/celia/Documents/VisualStudioCode/{dataset.upper()}/cellstatesappli/RNAmatrix.tsv')
    data = pd.read_csv(data_path, delimiter='\t', index_col=0, nrows=1)
    return data.columns.tolist()

def get_cells_barcodes(node, sparse_matrices):
    dataset = node.split('-')[-1].lower()
    nodes = node_lists[dataset]
    cells = cell_lists[dataset]
    node_index = nodes.index(node)
    matching_indices = sparse_matrices[dataset].row == node_index
    columns = sparse_matrices[dataset].col[matching_indices]
    corresponding_cells = [cells[col_index] for col_index in columns]
    return corresponding_cells

def process_list(lst):
    processed_mapping = {}
    for elem in lst:
        # Split and join parts except the first
        processed_part = '_'.join(elem.split('_')[1:])
        if processed_part not in processed_mapping:
            processed_mapping[processed_part] = []
        processed_mapping[processed_part].append(elem)
    return processed_mapping

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process hierarchical clustering and generate a pairwise matrix comparing clusters.')
    parser.add_argument('clustering_output_path', type=str, help='Path to clustering output file.')
    parser.add_argument('output_path', type=str, help='Path to save the output pairwise matrix CSV file.')
    parser.add_argument('sparse_matrix_dir', type=str, help='Directory containing the sparse matrix files')
    args = parser.parse_args()

    sparse_matrix_dir = Path(args.sparse_matrix_dir)
    data = ['peng', 'steele', 'lin', 'hwang', 'zhang', 'werba']
    
    # Load sparse matrices
    sparse_matrices = {}
    for dataset in data:
        sparse_matrix_path = sparse_matrix_dir / f'cell_mtx_{dataset}.npz'
        if sparse_matrix_path.exists():
            sparse_matrices[dataset] = sp.load_npz(sparse_matrix_path)
        else:
            print(f"Warning: Sparse matrix file not found for dataset {dataset} at {sparse_matrix_path}")

    # Initialize node lists
    node_lists = {}
    for dataset in data:
        node_lists[dataset] = create_dataset_node_list(dataset)

    # Initialize cell lists
    cell_lists = {}
    for dataset in data:
        cell_lists[dataset] = create_dataset_cell_list(dataset)

    # Load clustering output
    nodes_description = pd.read_csv(args.clustering_output_path, sep='\t')

    # Create cluster to node mapping
    cluster_node_mapping = {}
    for _, row in tqdm(nodes_description.iterrows(), total=nodes_description.shape[0], desc="Processing nodes"):
        node = row['Node']
        cluster = row['Cluster']
        if cluster not in cluster_node_mapping:
            cluster_node_mapping[cluster] = []
        cluster_node_mapping[cluster].append(node)

    # Create pairwise matrix
    cluster_names = list(cluster_node_mapping.keys())
    pairwise_matrix = pd.DataFrame(0, index=cluster_names, columns=cluster_names, dtype=float)

    # Populate pairwise matrix with comparisons
    for label in tqdm(cluster_names, desc="Processing clusters"):
        nodes = cluster_node_mapping.get(label, [])
        processed_lbl_nodes = process_list(nodes)

        for other_label in tqdm(cluster_names, desc=f"Comparing with {label}", leave=False):
            if label != other_label:
                other_nodes = cluster_node_mapping.get(other_label, [])
                processed_other_nodes = process_list(other_nodes)

                common_parts = set(processed_lbl_nodes.keys()).intersection(processed_other_nodes.keys())

                right_orientation = 0
                all_comparisons = 0

                for part in common_parts:
                    lbl_nodes = processed_lbl_nodes[part]
                    parent_nodes = processed_other_nodes[part]
                    for lbl in lbl_nodes:
                        lbl_cell_num = get_cells_barcodes(lbl, sparse_matrices)
                        for parent in parent_nodes:
                            parent_cell_num = get_cells_barcodes(parent, sparse_matrices)
                            if len(set(lbl_cell_num).intersection(parent_cell_num)) > 0:
                                all_comparisons += 1
                                if len(lbl_cell_num) < len(parent_cell_num):
                                    right_orientation += 1

                pairwise_matrix.loc[label, other_label] = right_orientation / all_comparisons if all_comparisons > 0 else 0

    # Save the pairwise matrix to a CSV file
    pairwise_matrix.to_csv(args.output_path)



# annotation_path = '/media/celia/data/colors.csv'
# condensed_matrix_path = '/media/celia/data/condensed_matrices/condensed_matrix_zhang_RBH02_cut08.csv'
# labels_path = '/media/celia/data/condensed_matrices/labels_zhang_RBH02_cut08.csv'
# pairwise_matrix.to_csv('/media/celia/data/condensed_matrices/clusters_comparisons3.csv')