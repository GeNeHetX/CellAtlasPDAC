import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path
from ete3 import Tree 
import scipy.sparse as sp 
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from tqdm import tqdm
import pickle
from scipy.sparse import csr_matrix



clustering_output =  pd.read_csv('/media/celia/data/atlas2/clustering_output_filtered.csv', sep = ',')

def extract_nodes(nodes_list, pattern_list):
    extracted_nodes = [node for node in nodes_list if any(pattern in node for pattern in pattern_list)]
    return extracted_nodes

datasets_name = ['peng' , 'steele', 'lin' , 'hwang', 'zhang', 'werba']
cell_mtx_dir = Path('/media/celia/data/cell_mtx/')

node_lists_ = {}
cells_lists_ = {}
sparse_matrices_ = {}

for dataset in datasets_name:
    sparse_data_path = cell_mtx_dir / f'cell_mtx_{dataset}_data.npz'
    if sparse_data_path.exists():
        data = np.load(sparse_data_path, allow_pickle=True)
        sparse_matrices_[dataset] = csr_matrix(data['sparse_matrix'].item())
        node_lists_[dataset] = data['node_index_map']
        cells_lists_[dataset] = data['cell_index_map']
    else:
        print(f"Warning: Sparse matrix file not found for dataset {dataset} at {sparse_data_path}")

def get_cell_barcodes(node):
    dataset = node.split('-')[-1]
    node_index_map = node_lists_[dataset]
    sparse_mtx = sparse_matrices_[dataset]
    cell_index_map = cells_lists_[dataset]
    if node not in node_index_map:
        raise ValueError(f"Node '{node}' not found in the index map")

    row_idx = np.where(node_index_map == node)[0][0]
    cells_indices = sparse_mtx[row_idx].nonzero()[1]
    cell_barcodes = cell_index_map[cells_indices]
    return cell_barcodes


cluster_node_mapping = {}
for _, row in clustering_output.iterrows():
    node = row['Node']
    cluster = row['Cluster']
    if cluster not in cluster_node_mapping:
        cluster_node_mapping[cluster] = []
    cluster_node_mapping[cluster].append(node)

# Precompute samples for each cluster and map samples to nodes for efficient filtering
cluster_samples = {
    clust: {node.split('_')[-1] for node in nodes}
    for clust, nodes in cluster_node_mapping.items()
}

# Create a sample-to-nodes mapping to speed up node filtering in `check_direct_link`
sample_to_nodes = defaultdict(list)
for clust, nodes in cluster_node_mapping.items():
    for node in nodes:
        sample = node.split('_')[-1]
        sample_to_nodes[sample].append(node)

# Cache cell barcodes for each node to avoid redundant calls to `get_cell_barcodes`
# node_cell_cache = {}
# def preload_node_cells():
#     def cache_node_cells(node):
#         node_cell_cache[node] = get_cell_barcodes(node)
    
#     all_nodes = [node for nodes in node_lists_.values() for node in nodes]
        
#     with ThreadPoolExecutor() as executor:
#         list(tqdm(executor.map(cache_node_cells, all_nodes), total = len(all_nodes), desc= 'Preloading nodes'))
    
#     with open('node_cell_cache.pkl', 'wb') as f:
#         pickle.dump(node_cell_cache, f)

# preload_node_cells()

def load_node_cell_cache():
    with open('node_cell_cache.pkl', 'rb') as f:
        return pickle.load(f)

node_cell_cache = load_node_cell_cache()


def get_cached_cells(node):
    if node not in node_cell_cache:
        node_cell_cache[node] = get_cell_barcodes(node)
    return node_cell_cache[node]

def check_direct_link(node1, node2, clust1, clust2, node1_cells, node2_cells):

    def get_nodes_to_check(sample, clust1, clust2):
        clusters_to_exclude = {clust1, clust2}

        nodes_to_check = [
        node for key,items in cluster_node_mapping.items()
        if key not in clusters_to_exclude
        for node in items if sample in node
        ]

        return nodes_to_check

    node1_nb_cells = len(node1_cells)
    node2_nb_cells = len(node2_cells)
    sample = node1.split("_")[1]
    nodes_to_check = get_nodes_to_check(sample, clust1, clust2)
    
    # Check potential intermediary nodes
    for node in nodes_to_check:
        node_cells = get_cached_cells(node)
        node_nb_cells = len(node_cells)
        if set(node1_cells).intersection(node_cells) and set(node2_cells).intersection(node_cells):
            if (min(node1_nb_cells, node2_nb_cells) <= node_nb_cells <= max(node1_nb_cells, node2_nb_cells)):
                return False
    return True

def calculate_branch_and_orientation(common_sample_nodes1, common_sample_nodes2, clust1, clust2):
    print(clust1, clust2)
    # Initialize counters
    same_branch = 0
    diff_branch = 0
    clust1to2 = 0
    clust2to1 = 0

    # Precompute cell barcodes for nodes in this sample
    node1_cells_list = [get_cached_cells(node1) for node1 in common_sample_nodes1]
    node2_cells_list = [get_cached_cells(node2) for node2 in common_sample_nodes2]

    # Iterate over nodes and calculate branch and orientation
    for node1, node1_cells in zip(common_sample_nodes1, node1_cells_list):
        for node2, node2_cells in zip(common_sample_nodes2, node2_cells_list):
            if set(node1_cells).intersection(node2_cells):
                same_branch += 1
                if check_direct_link(node1, node2, clust1, clust2, node1_cells, node2_cells):
                    if len(node1_cells) < len(node2_cells):
                        clust1to2 += 1
                    elif len(node1_cells) > len(node2_cells):
                        clust2to1 += 1
            else:
                diff_branch += 1

    return same_branch, diff_branch, clust1to2, clust2to1


def calculate_for_cluster_pair(clust1, clust2):
    result = {
        'same_sample': 0,
        'samples_clst1': 0,
        'samples_clst2': 0,
        'same_branch': 0,
        'diff_branch': 0,
        'orientation1to2': 0,
        'orientation2to1': 0
    }
    
    samples1 = cluster_samples[clust1]
    samples2 = cluster_samples[clust2]
    result['samples_clst1'] = len(samples1)
    result['samples_clst2'] = len(samples2)
    common_samples = samples1 & samples2
    result['same_sample'] = len(common_samples)

    for sample in common_samples:
        common_sample_nodes1 = extract_nodes(cluster_node_mapping[clust1], [sample])
        common_sample_nodes2 = extract_nodes(cluster_node_mapping[clust2], [sample])
        same_branch, diff_branch, clust1to2, clust2to1 = calculate_branch_and_orientation(
            common_sample_nodes1, common_sample_nodes2, clust1, clust2
        )
        result['same_branch'] += same_branch
        result['diff_branch'] += diff_branch
        result['orientation1to2'] += clust1to2
        result['orientation2to1'] += clust2to1

    return (clust1, clust2), result

def process_cluster_pair(pair):
    clust1, clust2 = pair
    return calculate_for_cluster_pair(clust1, clust2)

def parallel_process_clusters():
    cluster_pairs = [
        (clust1, clust2)
        for i, clust1 in enumerate(cluster_node_mapping.keys())
        for clust2 in list(cluster_node_mapping.keys())[i+1:]
    ]
    
    results = []
    with ProcessPoolExecutor(max_workers=18) as executor:
        results = list(executor.map(process_cluster_pair, cluster_pairs))

    # Aggregate results
    combined_data = []
    for (clust1, clust2), result in results:
        combined_data.append({
            'Cluster1': clust1,
            'Cluster2': clust2,
            'Cluster1 Samples': result['samples_clst1'],
            'Cluster2 Samples': result['samples_clst2'],
            'Number of Same Samples': result['same_sample'],
            'Same Branch Total': result['same_branch'],
            'Different Branch Total': result['diff_branch'],
            'Cluster2 Parent': result['orientation1to2'],
            'Cluster1 Parent': result['orientation2to1'],
        })
    
    return combined_data

# Run the parallel processing and save results
combined_data = parallel_process_clusters()
df = pd.DataFrame(combined_data)
df.to_csv('/media/celia/data/atlas2/cluster_comparison_direct_links.csv', index=False)
