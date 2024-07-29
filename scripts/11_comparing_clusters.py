import pandas as pd
import numpy as np
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette
import argparse
from ete3 import Tree
import scipy.sparse as sp

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

def initialize_data(condensed_matrix_path, labels_path):

    condensed_matrix = np.loadtxt(condensed_matrix_path, delimiter=',')
    with open(labels_path, 'r') as f:
        labels = [line.strip() for line in f]

    Z = linkage(condensed_matrix, 'average')
    set_link_color_palette([f"C{i}" for i in range(1, 100000)])
    d = dendrogram(Z, labels=labels, color_threshold=0.8, above_threshold_color='black', no_plot=True)

    trtm_status_hwang = {'U': 'untreated', 'T': 'treated'}
    trtm_status = {'peng': 'untreated', 'steele': 'untreated', 'lin': 'untreated', 'zhang': 'untreated'}

    nodes = []
    clusters = []
    treatments = []
    datasets = []

    for i, node in enumerate(d['ivl']):
        nodes.append(node)
        clusters.append(d['leaves_color_list'][i])
        if node.split('-')[-1] == 'hwang':
            trtm = trtm_status_hwang.get(node.split('_')[2][0], 'black')
        else:
            trtm = trtm_status.get(node.split('-')[-1])
        treatments.append(trtm)
        datasets.append(node.split('-')[-1])

    nodes_description = pd.DataFrame({
        'Node': nodes,
        'Cluster': clusters,
        'Dataset': datasets,
        'Treatment': treatments
    })

    return nodes_description

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process hierarchical clustering and generate a pairwise matrix comparing clusters.')
    parser.add_argument('condensed_matrix_path', type=str, help='Path to the condensed matrix CSV file.')
    parser.add_argument('labels_path', type=str, help='Path to the labels CSV file.')
    parser.add_argument('output_path', type=str, help='Path to save the output pairwise matrix CSV file.')
    parser.add_argument("--sparse_matrix_dir", type=str, required=True, help="Directory containing the sparse matrix files")
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

    nodes_description = initialize_data(args.condensed_matrix_path, args.labels_path)

    cluster_cell_mapping = {}

    for i, (_, row) in enumerate(nodes_description.iterrows()):
        node = row['Node']
        cluster = row['Cluster']
        print(f'1/2 \t {i}/{len(nodes_description.index)}')
        cell_barcodes = get_cells_barcodes(node, sparse_matrices)
        
        if cluster not in cluster_cell_mapping:
            cluster_cell_mapping[cluster] = []
        cluster_cell_mapping[cluster].extend(cell_barcodes)

    cluster_names = list(cluster_cell_mapping.keys())
    pairwise_matrix = pd.DataFrame(0, index=cluster_names, columns=cluster_names, dtype=float)


    for i, (cluster, cell_barcodes) in enumerate(cluster_cell_mapping.items()):
        print(f'2/2 {i}/{len(cluster_cell_mapping.keys())}')
        for other_cluster in cluster_names:
            if other_cluster != cluster:
                cells_in_other_cluster = sum(1 for cell in cell_barcodes if cell in cluster_cell_mapping[other_cluster])
                pairwise_matrix.loc[cluster, other_cluster] = cells_in_other_cluster / len(cell_barcodes) if len(cell_barcodes) > 0 else 0


    pairwise_matrix.to_csv(args.output_path)



# annotation_path = '/media/celia/data/colors.csv'
# condensed_matrix_path = '/media/celia/data/condensed_matrices/condensed_matrix_zhang_RBH02_cut08.csv'
# labels_path = '/media/celia/data/condensed_matrices/labels_zhang_RBH02_cut08.csv'
# pairwise_matrix.to_csv('/media/celia/data/condensed_matrices/clusters_comparisons3.csv')