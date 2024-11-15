import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette
import os
from tqdm.auto import tqdm
import argparse

def read_and_prepare_data(input_path, distance_column):
    """
    Read input data and prepare the distance column.

    Parameters:
    input_path (str): Path to the input CSV file.
    distance_column (str): Column name containing the distance metric.

    Returns:
    data (pd.DataFrame): DataFrame with prepared distance column.
    nodes (list): List of unique nodes.
    node_indices (dict): Dictionary mapping nodes to indices.
    """
    print('Preparing data')
    data = pd.read_csv(input_path, sep='\t')
    data['Distance'] = 1 - data[distance_column]
    nodes = list(set(data['Node1']).union(data['Node2']))
    node_indices = {node: i for i, node in enumerate(nodes)}
    return data, nodes, node_indices


def create_distance_matrix(data, nodes, node_indices):
    """
    Create a distance matrix from the data.

    Parameters:
    data (pd.DataFrame): Input data.
    nodes (list): List of unique nodes.
    node_indices (dict): Dictionary mapping nodes to indices.

    Returns:
    distance_matrix (np.array): Distance matrix.
    """
    n = len(nodes)
    print(f"Number of nodes: {n}")
    print('Creating distance matrix')
    
    # Initialize the distance matrix with ones and set the diagonal to zero
    distance_matrix = np.ones((n, n))
    np.fill_diagonal(distance_matrix, 0)
    
    # Map 'Node1' and 'Node2' to their corresponding indices in the distance matrix
    indices_1 = data['Node1'].map(node_indices)
    indices_2 = data['Node2'].map(node_indices)
    
    # Assign the distances 
    distance_matrix[indices_1, indices_2] = data['Distance'].values
    distance_matrix[indices_2, indices_1] = data['Distance'].values
    
    return distance_matrix

def perform_clustering(output_dir, distance_matrix, nodes, linkage_method, color_threshold, min_cluster_size):
    """
    Perform hierarchical clustering and filter small clusters.

    Parameters:
    distance_matrix (np.array): Distance matrix.
    nodes (list): List of unique nodes.
    linkage_method (str): Linkage method for hierarchical clustering.
    color_threshold (float): Threshold for coloring clusters in the dendrogram.
    min_cluster_size (int): Minimum size of clusters to keep.

    Returns:
    filtered_nodes (list): List of nodes in valid clusters.
    """
    print('Condensing distance matrix')
    condensed_matrix = squareform(distance_matrix)

    print('Computing linkage')
    Z = linkage(condensed_matrix, linkage_method)

    print('Computing dendrogram')
    set_link_color_palette([f"C{i}" for i in range(1, 100000)])
    dn = dendrogram(Z, labels=nodes, color_threshold=color_threshold, no_plot=True)

    clusters = pd.DataFrame(data={'Node': dn['ivl'], 'Cluster': dn['leaves_color_list']})
    filtered_clusters = clusters[clusters['Cluster'] != 'C0']
    cluster_counts = filtered_clusters['Cluster'].value_counts()
    valid_clusters = cluster_counts[cluster_counts >= min_cluster_size].index
    filtered_clusters = filtered_clusters[filtered_clusters['Cluster'].isin(valid_clusters)]
    
    nodes = []
    clusters = []
    for i, node in enumerate(dn['ivl']):
        nodes.append(node)
        clusters.append(dn['leaves_color_list'][i])

    res = pd.DataFrame({'Node': nodes, 'Cluster': clusters})

    # res.to_csv(os.path.join(output_dir, 'first_clustering_output.csv'), index=False, sep='\t')

    # save_results(output_dir, condensed_matrix, nodes)

    return res

def save_results(output_dir, filtered_condensed_matrix, nodes):
    """
    Save the condensed matrix and node labels to files.

    Parameters:
    output_dir (str): Directory to save the output files.
    filtered_condensed_matrix (np.array): Condensed distance matrix.
    nodes (list): List of nodes.
    """
    os.makedirs(output_dir, exist_ok=True)
    condensed_matrix_path = os.path.join(output_dir, 'condensed_mtx.csv')
    np.savetxt(condensed_matrix_path, filtered_condensed_matrix, delimiter=',')

    labels_path = os.path.join(output_dir, 'labels.csv')
    with open(labels_path, 'w') as f:
        for label in nodes:
            f.write(f"{label}\n")

def main(input_path, output_dir, distance_column='Jaccard Score', linkage_method='average', cut_tree_threshold=0.8, min_cluster_size=3):
    """
    Main function to perform hierarchical clustering on node data and save filtered results.

    Parameters:
    input_path (str): Path to the input CSV file.
    output_dir (str): Directory to save the output files.
    distance_column (str): Column name containing the distance metric.
    linkage_method (str): Linkage method for hierarchical clustering.
    cut_tree_threshold (float): Threshold for clusters in the dendrogram.
    min_cluster_size (int): Minimum size of clusters to keep.
    """
    # Step 1: Read and prepare data
    try:
        data, nodes, node_indices = read_and_prepare_data(input_path, distance_column)
    except Exception as e:
        print(f"Error reading or preparing data: {e}")
        return
    
    # Step 2: Create a distance matrix
    try:
        distance_matrix = create_distance_matrix(data, nodes, node_indices)
        print("Distance Matrix:\n", distance_matrix)
    except Exception as e:
        print(f"Error creating distance matrix: {e}")
        return

    # Step 3: Perform hierarchical clustering
    try:
        clustering_result = perform_clustering(
            output_dir=output_dir,
            distance_matrix=distance_matrix,
            nodes=nodes,
            linkage_method=linkage_method,
            color_threshold=cut_tree_threshold,
            min_cluster_size=min_cluster_size
        )
    except Exception as e:
        print(f"Error during clustering: {e}")
        return
    
    # Step 4: Filter clusters
    if clustering_result is not None:
        filtered_clusters = clustering_result[clustering_result['Cluster'] != 'C0']

        # Count the number of nodes in each cluster
        cluster_counts = filtered_clusters['Cluster'].value_counts()

        # Filter clusters by minimum size
        valid_clusters = cluster_counts[cluster_counts >= min_cluster_size].index
        filtered_clusters = filtered_clusters[filtered_clusters['Cluster'].isin(valid_clusters)]

        # Rename clusters sequentially (e.g., C1, C2, C3, ...)
        cluster_mapping = {old_cluster: f'C{i+1}' for i, old_cluster in enumerate(sorted(filtered_clusters['Cluster'].unique()))}
        filtered_clusters['Cluster'] = filtered_clusters['Cluster'].map(cluster_mapping)
   
        # Step 5: Save the filtered clustering results
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, 'clustering_output.csv')
        filtered_clusters.to_csv(output_file, index=False, sep='\t')
        
        print(f"Filtered clustering results saved to: {output_file}")
    else:
        print("No valid clusters found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform hierarchical clustering on node data and save filtered results.')
    parser.add_argument('input_path', type=str, help='Path to the input CSV file : filtered computed similarities.')
    parser.add_argument('output_dir', type=str, help='Directory to save the output files: condensed matrix with its labels.')
    parser.add_argument('--distance_column', type=str, default='Jaccard Score', help='Column name containing the distance metric.')
    parser.add_argument('--linkage_method', type=str, default='average', help='Linkage method for hierarchical clustering.')
    parser.add_argument('--cut_tree_threshold', type=float, default=0.8, help='Threshold for clusters in the dendrogram.')
    parser.add_argument('--min_cluster_size', type=int, default=3, help='Minimum size of clusters to keep.')

    args = parser.parse_args()

    main(args.input_path, args.output_dir, args.distance_column, args.linkage_method, args.cut_tree_threshold, args.min_cluster_size)

# python script_name.py /path/to/input/full_data.csv /path/to/output_dir --distance_column 'Jaccard Score' --linkage_method 'average' --cut_tree_threshold 0.8 --min_cluster_size 3
# input is similarity file filtered