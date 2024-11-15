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

def has_enough_cells(node, dataset, row_sums):
    dataset = dataset.lower()
    node_index = np.where(node_lists[dataset] == node)[0][0]
    row_sum = row_sums[node_index]
    return row_sum >= 10


def process_genes(node, gene_mtx_subset_path, dataset, subset, num_cells):
    if has_enough_cells(node, dataset, num_cells[dataset]):
        gene_mtx_data = np.load(gene_mtx_subset_path, allow_pickle=True)

        try:
            sparse_matrix_over = gene_mtx_data['sparse_matrix_over'].item()
            sparse_matrix_under = gene_mtx_data['sparse_matrix_under'].item()
            gene_index_map_over = dict(gene_mtx_data['gene_index_map_over'])
            gene_index_map_under = dict(gene_mtx_data['gene_index_map_under'])
            node_index_map = dict(gene_mtx_data['node_index_map'])
        except KeyError as e:
            raise KeyError(f"Missing key in data file: {e}")
        
        n = node.split('-')[0].split('_')[0]
        row_idx = node_index_map[n]

        over_gene_indices = sparse_matrix_over[row_idx].nonzero()[1]
        over_expressed_genes = np.array(list(gene_index_map_over.keys()))[over_gene_indices]

        under_gene_indices = sparse_matrix_under[row_idx].nonzero()[1]
        under_expressed_genes = np.array(list(gene_index_map_under.keys()))[under_gene_indices]

        expression = {}

        # Assign +2 for over-expressed genes
        for gene in over_expressed_genes:
            expression[gene] = +2
        
        # Assign -2 for under-expressed genes
        for gene in under_expressed_genes:
            expression[gene] = -2

        return node, expression
    return None, None

def create_dict(dataset, subset, num_workers, num_cells, gene_mtx_path):
    genes_dict = {}
    gene_mtx_subset_path = Path(f'{gene_mtx_path}/{dataset}/{dataset}_{subset}_data.npz')
    nodes = node_lists[dataset]
    nodes_sample = [node for node in nodes if node.split('_')[1] == f'{subset}-{dataset}']

    with Pool(num_workers) as pool:
        results = pool.starmap(process_genes, [(node, gene_mtx_subset_path, dataset, subset, num_cells) for node in nodes_sample])


    for node, expression in results:
        if node and expression:
            genes_dict[node] = expression

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

def main(output_path, cell_mtx_dir, gene_mtx_dir, samples_dir, num_processes):
    # Define all paths
    output_path = Path(output_path)
    cell_mtx_dir = Path(cell_mtx_dir)
    sample_dir = Path(samples_dir)
    datasets_name = ['peng' , 'steele', 'lin', 'hwang', 'zhang', 'werba']
    
    # Load sparse matrices and node lists
    global node_lists

    node_lists = {}
    sparse_matrices = {}

    for dataset in datasets_name:
        sparse_data_path = cell_mtx_dir / f'cell_mtx_{dataset}_data.npz'
        if sparse_data_path.exists():
            data = np.load(sparse_data_path, allow_pickle=True)
            sparse_matrices[dataset] = data['sparse_matrix'].item()
            node_lists[dataset] = data['node_index_map']
        else:
            print(f"Warning: Sparse matrix file not found for dataset {dataset} at {sparse_data_path}")

    # Compute row sums
    row_sums = {
        dataset: sparse_matrix.sum(axis=1).A.ravel()
        for dataset, sparse_matrix in sparse_matrices.items()
    }
    
    datasets = {}
    for dataset in datasets_name:
        file_path = sample_dir / f'{dataset}_samples.txt'
        
        if file_path.exists():
            with open(file_path, 'r') as f:
                samples = [line.strip() for line in f.readlines()]
                datasets[dataset] = samples
        else:
            print(f"Warning: File {file_path} not found!")

    # Create genes dictionary
    global get_genes_dict
    get_genes_dict = {}
    for dataset, subsets in datasets.items():
        dict_pbar = tqdm(total=len(subsets), desc=f"Processing subsets for {dataset}")
        genes_dict_subset = {}
        for subset in subsets:
            genes_dict = create_dict(dataset, subset, num_processes, row_sums, gene_mtx_dir)
            genes_dict_subset[subset] = genes_dict
            dict_pbar.update(1)
        get_genes_dict[dataset] = genes_dict_subset
        dict_pbar.close()

    pool = multiprocessing.Pool(processes=num_processes)

    # Create list of subset combinations
    prefixed_subsets = {f"{dataset}--{subset}": subset for dataset, subsets in datasets.items() for subset in subsets}
    all_subsets_with_prefix = list(prefixed_subsets.keys())
    subset_combinations = list(itertools.combinations(all_subsets_with_prefix, 2))

    sample_to_avoid = ['Case2.ZC'] # ADD samples to not take into account

    # Prepare tasks for multiprocessing
    tasks = []
    total_tasks = 0
    for comb1, comb2 in subset_combinations:
        dataset1 = comb1.split('--')[0]
        dataset2 = comb2.split('--')[0]
        subset1 = comb1.split('--')[1]
        subset2 = comb2.split('--')[1]

        if (subset1 not in sample_to_avoid) & (subset2 not in sample_to_avoid):
            if dataset1 == 'peng' and dataset2 == 'peng':
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
    parser.add_argument("--cell_mtx_dir", type=str, required=True, help="Directory containing the sparse cell matrix files")
    parser.add_argument("--gene_mtx_dir", type=str, required=True, help="Directory containing the sparse gene matrix files")
    parser.add_argument("--samples_dir", type=str, required=True, help="Directory containing files of sample id for each dataset")
    parser.add_argument("--num_processes", type=int, default=10, help="Number of processes for multiprocessing")
    args = parser.parse_args()

    main(output_path=args.output_path,
        cell_mtx_dir=args.cell_mtx_dir,
        gene_mtx_dir=args.gene_mtx_dir,
        samples_dir = args.samples_dir,
        num_processes=args.num_processes)


# python3 07_similarity_rbh.py --output_path /media/celia/data/nodes/ --cell_mtx_dir /media/celia/data/cell_mtx/  --gene_mtx_dir /media/celia/data/gene_mtx/ --samples_dir /home/celia/Documents/VisualStudioCode/samples/ --num_processes 2


