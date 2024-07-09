import argparse
import csv
import multiprocessing
import itertools
from pathlib import Path
import pandas as pd
import numpy as np
from ete3 import Tree
import os
import scipy.sparse as sp
from multiprocessing import Pool
from tqdm import tqdm
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

def create_node_ids(subset, dataset, cellstates_base_path):
    tree_path = f"{cellstates_base_path}{dataset.upper()}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"
    t = Tree(tree_path, format=1)
    return [f"{node.name}_{subset}-{dataset}" for node in t.traverse() if not node.is_root()]


def create_dataset_node_list(dataset, datasets, cellstates_base_path):
    dataset_node_list = []
    for subset in datasets[dataset]:
        dataset_node_list.extend(create_node_ids(subset, dataset, cellstates_base_path))
    return dataset_node_list

def has_enough_cells(node, dataset, row_sums):
    dataset = dataset.lower()
    node_index = node_lists[dataset].index(node)
    row_sum = row_sums[node_index]
    return row_sum >= 10

def process_file(file_path, dataset, subset, num_cells):
    node = file_path.stem  # Extract node name from file name
    if has_enough_cells(f'{node}_{subset}-{dataset}', dataset, num_cells[dataset]):
        try:
            marker_file = pd.read_csv(file_path, usecols=[0,2,5])
            marker_file.rename(columns={marker_file.columns[0]: 'genes'}, inplace=True)
            if dataset == 'hwang':
                marker_file['genes'] = marker_file['genes'].str.split('\t').str[0]
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            return None, None
        
        marker_file_filtered = marker_file[marker_file['p_val_adj'] <= 0.05]
        
        if marker_file_filtered.empty:
            return None, None

        genes_logfc_filtered = marker_file_filtered.set_index('genes')['avg_log2FC'].to_dict()
        new_node = f'{node}_{subset}-{dataset}'
        return new_node, genes_logfc_filtered
    return None, None

def create_dict(dataset, subset, num_workers, num_cells, markers_base_path):
    genes_dict = {}
    markers_folder = Path(f'{markers_base_path}/{dataset}/')
    subset_folder = markers_folder / subset
    files = list(subset_folder.glob('*.csv'))

    with Pool(num_workers) as pool:
        results = pool.starmap(process_file, [(file_path, dataset, subset, num_cells) for file_path in files])

    for node, genes_logfc_filtered in results:
        if node and genes_logfc_filtered:
            genes_dict[node] = genes_logfc_filtered

    return genes_dict

def compute_jaccard(dict1, dict2):
    nodes1 = sorted([node for node in dict1.keys()])
    nodes2 = sorted([node for node in dict2.keys()])
    jaccard_matrix = np.zeros((len(nodes1), len(nodes2)))

    for i, node1 in enumerate(nodes1):
        genes1 = dict1.get(node1, {})
        genes1_positive = {gene for gene, logfc in genes1.items() if logfc > 0}
        genes1_negative = {gene for gene, logfc in genes1.items() if logfc < 0}
        for j, node2 in enumerate(nodes2):
            genes2 = dict2.get(node2, {})
            genes2_positive = {gene for gene, logfc in genes2.items() if logfc > 0}
            genes2_negative = {gene for gene, logfc in genes2.items() if logfc < 0}
            intersection_positive = len(genes1_positive.intersection(genes2_positive))
            intersection_negative = len(genes1_negative.intersection(genes2_negative))
            union_genes = len(set(dict1[node1].keys()).union(set(dict2[node2].keys())))
            index = (intersection_positive + intersection_negative) / union_genes
            jaccard_matrix[i, j] = index

    jaccard_df = pd.DataFrame(jaccard_matrix, index=nodes1, columns=nodes2)
    return jaccard_df

def find_matching_pairs(df):
    matching_pairs = []
    for row_index, row in df.iterrows():
        max_col = row.idxmax()
        max_row_index = df[max_col].idxmax()
        max_value = row[max_col]
        if max_row_index == row_index:
            pair1 = (max_col, row_index, max_value)
            pair2 = (row_index, max_col, max_value)
            if pair1 not in matching_pairs and pair2 not in matching_pairs:
                matching_pairs.append(pair1)
    return matching_pairs

def get_output_file_path(subset1, dataset1, subset2, dataset2, base_path):
    return base_path / f'graph_nodes_data_{subset1}_{dataset1}_{subset2}_{dataset2}.csv'

def process_combination(subset1, dataset1, subset2, dataset2, output_file_path):
    print(f"Processing {subset1} from {dataset1} and {subset2} from {dataset2}")

    genes_dict1 = get_genes_dict.get(dataset1, {}).get(subset1)
    genes_dict2 = get_genes_dict.get(dataset2, {}).get(subset2)

    if genes_dict1 is None or genes_dict2 is None:
        print(f"Error: Genes dictionary not found for one or both subsets: {subset1}, {subset2}")
        return

    jaccard_matrix = compute_jaccard(genes_dict1, genes_dict2)

    rbh = find_matching_pairs(jaccard_matrix)
    del genes_dict1
    del genes_dict2
    with open(output_file_path, 'a', newline='') as csvfile:
        fieldnames = ['Node1', 'Node2', 'Jaccard Score', 'RBH', 'Same Dataset']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for node1 in jaccard_matrix.index:
            for node2 in jaccard_matrix.columns:
                jaccard_score = jaccard_matrix.loc[node1, node2]
                if jaccard_score != 0:
                    node1_dataset = node1.split('-')[1]
                    node2_dataset = node2.split('-')[1]
                    is_same_dataset = node1_dataset == node2_dataset
                    is_rbh = (node1, node2, jaccard_score) in rbh or (node2, node1, jaccard_score) in rbh
                    writer.writerow({'Node1': node1, 'Node2': node2, 'Jaccard Score': jaccard_score,
                                        'RBH': is_rbh, 'Same Dataset': is_same_dataset})
    print(f"Processing for {subset1} and {subset2} done!")

def main(output_path, cellstates_base_path, sparse_matrix_dir, markers_base_path, num_processes):
    # Define all paths
    output_path = Path(output_path)
    sparse_matrix_dir = Path(sparse_matrix_dir)
    datasets = ['peng', 'steele', 'lin', 'hwang', 'zhang']
    
    # Load sparse matrices
    sparse_matrices = {}
    for dataset in datasets:
        sparse_matrix_path = sparse_matrix_dir / f'cell_mtx_{dataset}.npz'
        if sparse_matrix_path.exists():
            sparse_matrices[dataset] = sp.load_npz(sparse_matrix_path)
        else:
            print(f"Warning: Sparse matrix file not found for dataset {dataset} at {sparse_matrix_path}")

    # Compute row sums
    row_sums = {
        dataset: sparse_matrix.sum(axis=1).A.ravel()
        for dataset, sparse_matrix in sparse_matrices.items()
    }


    # Define subsets for each dataset
    datasets = {
        'peng': ['T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22'],
        'steele': ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11A", "11B", "12", "13", "15", "16"],
        'lin': ["P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "MET01", "MET02", "MET03", "MET04", "MET05", "MET06"],
        'hwang': ['003_U15', '004_U18', '007_U2', '008_T22', '009_T20', '010T_U12', '011_U9', '2083_T4', '2100_T2', '2223_T17',
                '2229_T18', '2276_U4', '2364_U5', '2376_U1', '2443_U10', '2462_T16', '2490_U11', '2498_U7', '2507_T25', '2523_U3',
                    '2540_T6', '2591_U13', '2603_U6', '2626_U16', '2634_T5', '2664_U14', '2667_T23', '2668_T7', '2675_T13', 'MGH2076_T9',
                    'MGH2101_T8', 'MGH2381_T19', 'MGH2675_T13', 'MGHR1_T14', 'MGHR2_T15', 'MGHR3_T3', 'MGHR5_T1', 'MGHR6_T10', 'MGHR7_T12', 'MGHR8_T24',
                        'MGHR9_T21', 'MGHR11_T11', 'MGHR16_U17', 'MGHR17_U8'],
        'zhang': ['Case1.YF', 'Case1.ZY', 'Case2.YF', 'Case2.ZC', 'Case2.ZY', 'Case3.YF', 'Case3.ZY', 'Case4.ZY']
    }

    # Create node lists
    global node_lists
    node_lists = {}
    for dataset in datasets:
        node_lists[dataset] = create_dataset_node_list(dataset, datasets, cellstates_base_path)

    # Create genes dictionary
    global get_genes_dict
    get_genes_dict = {}
    for dataset, subsets in datasets.items():
        dict_pbar = tqdm(total=len(subsets), desc=f"Processing subsets for {dataset}")
        genes_dict_subset = {}
        for subset in subsets:
            genes_dict = create_dict(dataset, subset, 20, row_sums, markers_base_path)
            genes_dict_subset[subset] = genes_dict
            dict_pbar.update(1)
        get_genes_dict[dataset] = genes_dict_subset
        dict_pbar.close()

    pool = multiprocessing.Pool(processes=num_processes)

    # Create list of subset combinations
    prefixed_subsets = {f"{dataset}--{subset}": subset for dataset, subsets in datasets.items() for subset in subsets}
    all_subsets_with_prefix = list(prefixed_subsets.keys())
    subset_combinations = list(itertools.combinations(all_subsets_with_prefix, 2))

    # Prepare tasks for multiprocessing
    tasks = []
    total_tasks = 0
    for comb1, comb2 in subset_combinations:
        dataset1 = comb1.split('--')[0]
        dataset2 = comb2.split('--')[0]
        subset1 = comb1.split('--')[1]
        subset2 = comb2.split('--')[1]
        output_file_path = get_output_file_path(subset1, dataset1, subset2, dataset2, output_path)
        if not os.path.exists(output_file_path):
            tasks.append((subset1, dataset1, subset2, dataset2, output_file_path))
            total_tasks += 1
        else:
            print(f"Output file '{output_file_path}' already exists. Skipping processing.")
    
    pool.starmap_async(process_combination, tasks)
    pool.close()
    pool.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process computing similarity")
    parser.add_argument("--output_path", type=str, required=True, help="Folder path for output files")
    parser.add_argument("--cellstates_base_path", type=str, required=True, help="Base path for all cellstates outputs (including newick tree files)")
    parser.add_argument("--sparse_matrix_dir", type=str, required=True, help="Directory containing the sparse matrix files")
    parser.add_argument("--markers_base_path", type=str, required=True, help="Base path for marker files")
    parser.add_argument("--num_processes", type=int, default=10, help="Number of processes for multiprocessing")
    args = parser.parse_args()

    main(output_path=args.output_path, cellstates_base_path=args.cellstates_base_path, sparse_matrix_dir=args.sparse_matrix_dir,
         markers_base_path=args.markers_base_path, num_processes=args.num_processes)


# python3 07_similarity_RBH.py --output_path /media/celia/data/nodes/     --sparse_matrix_dir /media/celia/data/     
# --cellstates_base_path /home/celia/Documents/VisualStudioCode/     --markers_base_path /home/celia/Documents/annotation/findmarkers/     
# --num_processes 2


